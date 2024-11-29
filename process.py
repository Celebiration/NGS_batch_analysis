import pandas as pd
import numpy as np
import os
import sys
import subprocess
import re
from Bio import SeqIO
from collections import defaultdict
from collections import Counter
lines='-----------------------------'
print(pd.__version__)
tmpbase=str(sys.argv[2])
outbase='/var/www/platform/public/tmp/ngs/result/'+str(sys.argv[1])
if tmpbase[-1] != '/':
	tmpbase=tmpbase+'/'
if outbase[-1] != '/':
	outbase=outbase+'/'
if os.path.exists(tmpbase+'mapping.tsv'):
	mapping=pd.read_csv(tmpbase+'mapping.tsv',sep='\t')
elif os.path.exists(tmpbase+'mapping.txt'):
	mapping=pd.read_csv(tmpbase+'mapping.txt',sep='\t')
elif os.path.exists(tmpbase+'mapping.csv'):
	mapping=pd.read_csv(tmpbase+'mapping.csv')
else:
	print(os.listdir(tmpbase))
	# subprocess.run(args=["python","/var/www/platform/app/scripts/ngs/taskfailed.py",sys.argv[1],sys.argv[3]])
	raise ValueError('\n# <font color=#FF0000>❌mapping.txt, mapping.tsv or mapping.csv not found.</font>\n')
mapping['job']=list(map(lambda x:str(x).lower(),mapping['job']))
mapping['sample_id']=mapping['sample_id'].astype(str)
header=mapping[[False]*len(mapping)]
skipped=[]

#将sample_id和control中的特殊字符替换为'_'
mapping['sample_id']=[re.sub(r'[^a-zA-Z0-9_-]', '_', i) for i in mapping['sample_id']]
mapping['control']=[(lambda x: re.sub(r'[^a-zA-Z0-9_-]', '_', x) if isinstance(x,str) and len(x)>0 else x)(i) for i in mapping['control']]

#去除adapter中的空格
mapping['adapter_5']=[i.replace(' ','') if isinstance(i,str) and len(i)>0 else i for i in mapping['adapter_5']]
mapping['adapter_3']=[i.replace(' ','') if isinstance(i,str) and len(i)>0 else i for i in mapping['adapter_3']]

#检查是否有重复的样本名
if any(mapping[['sample_id','job']].duplicated()):
	raise ValueError('\n# <font color=#FF0000>❌同一job下有重复的sample_id，请检查。</font>\n')


#首先将混合样本进行demultiplex
mixed=mapping[mapping['mixed']=='yes']
unmixed=mapping[mapping['mixed']!='yes']
if len(mixed)>0:
	print('\n# '+'开始拆分混合样本'+'\n')
	new_mixed=header
	mixed_files=np.unique(mixed['file_name'])
	for i,files in enumerate(mixed_files):
		if files.count(';')==0:#为单端测序
			data=mapping[mapping['file_name']==files]
			f=open(tmpbase+'demult_'+str(i+1)+'_F.fa','w')
			r=open(tmpbase+'demult_'+str(i+1)+'_R.fa','w')
			for ii in range(len(data)):
				f.write('>'+list(data['sample_id'])[ii]+'_f'+'\n'+list(data['adapter_5'])[ii][:15]+'\n')#barcode长度不超过15
				r.write('>'+list(data['sample_id'])[ii]+'_r'+'\n'+list(data['adapter_3'])[ii][:15]+'\n')
			f.close()
			r.close()
			R1=files.strip('.gz')
			print('\n## '+'开始对文件：%s进行拆分：%s' %(R1, str(list(data['sample_id'])))+'\n')
			val=os.system('cutadapt -j 8 --revcomp -g ^file:%s -a ^file:%s --pair-adapters -e 0.1 -o "%s{name}.fq" --action none --untrimmed-output %s %s' %(tmpbase+'demult_'+str(i+1)+'_F.fa',tmpbase+'demult_'+str(i+1)+'_R.fa',tmpbase+'fq/', tmpbase+'fq/'+R1+'.untrimmed.fq',tmpbase+'fq/'+R1))
			if val!=0:
				raise ValueError('\n### <font color=#FF0000>❌demultiplex时，cutadapt指令无法成功运行。</font>\n')
			print('完成。')
			#如果拆分出的reads比例<80%，则警告：
			before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+R1),shell=True).decode('utf-8').strip())
			after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+R1+'.untrimmed.fq'),shell=True).decode('utf-8').strip())
			if after/before > 0.2:
				print("\n## <font color=#D89F12>⚠️警告：未拆分出的reads比例为{:.2%}</font>\n".format(after/before))

			valid=[]
			for j,sample in enumerate(data['sample_id']):
				val=os.system('cat %s_f.fq %s_r.fq > %s_combined.fq' %(tmpbase+'fq/'+sample,tmpbase+'fq/'+sample,tmpbase+'fq/'+sample))
				#如果文件不存在，则删除条目并加入skipped
				if val!=0:
					print('\n## <font color=#D89F12>⚠️警告：跳过%s，因为混合样本%s中没有其对应的reads。</font>\n' %(sample,files))
					skipped.append(sample)
					continue
				_=os.system('rm %s_f.fq %s_r.fq' %(tmpbase+'fq/'+sample,tmpbase+'fq/'+sample))
				valid.append(j)
			data['file_name']=data['sample_id']+'_combined.fq'
			data=data.iloc[valid,]
			new_mixed=pd.concat([new_mixed,data])
		elif files.count(';')==1:#为双端测序
			data=mapping[mapping['file_name']==files]
			f=open(tmpbase+'demult_'+str(i+1)+'_F.fa','w')
			r=open(tmpbase+'demult_'+str(i+1)+'_R.fa','w')
			for ii in range(len(data)):
				f.write('>'+list(data['sample_id'])[ii]+'_f'+'\n'+list(data['adapter_5'])[ii][:15]+'\n'+'>'+list(data['sample_id'])[ii]+'_r'+'\n'+list(data['adapter_3'])[ii][:15]+'\n')#barcode长度不超过15
				r.write('>'+list(data['sample_id'])[ii]+'_r'+'\n'+list(data['adapter_3'])[ii][:15]+'\n'+'>'+list(data['sample_id'])[ii]+'_f'+'\n'+list(data['adapter_5'])[ii][:15]+'\n')
			f.close()
			r.close()
			R1,R2=files.split(';')
			R1=R1.strip('.gz')
			R2=R2.strip('.gz')
			print('\n## '+'开始对文件：%s进行拆分：%s' %(files, str(list(data['sample_id'])))+'\n')
			val=os.system('cutadapt -j 8 -g ^file:%s -G ^file:%s --pair-adapters -e 0.1 -o "%s{name}.1.fq" -p "%s{name}.2.fq" --action none --untrimmed-output %s --untrimmed-paired-output %s %s %s' %(tmpbase+'demult_'+str(i+1)+'_F.fa',tmpbase+'demult_'+str(i+1)+'_R.fa',tmpbase+'fq/',tmpbase+'fq/',tmpbase+'fq/'+R1+'_'+R2+'.untrimmed_1.fq',tmpbase+'fq/'+R1+'_'+R2+'.untrimmed_2.fq',tmpbase+'fq/'+R1,tmpbase+'fq/'+R2))
			if val!=0:
				raise ValueError('\n### <font color=#FF0000>❌demultiplex时，cutadapt指令无法成功运行。</font>\n')
			print('完成。')
			#如果拆分出的reads比例<80%，则警告：
			before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+R1),shell=True).decode('utf-8').strip())
			after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+R1+'_'+R2+'.untrimmed_1.fq'),shell=True).decode('utf-8').strip())
			if after/before > 0.2:
				print("\n## <font color=#D89F12>⚠️警告：未拆分出的reads比例为{:.2%}</font>\n".format(after/before))
			
			valid=[]
			for j,sample in enumerate(data['sample_id']):
				val=os.system('cat %s_f.1.fq %s_r.2.fq > %s_combined.1.fq' %(tmpbase+'fq/'+sample,tmpbase+'fq/'+sample,tmpbase+'fq/'+sample))
				#如果文件不存在，则删除条目并加入skipped
				if val!=0:
					print('\n## <font color=#D89F12>⚠️警告：跳过%s，因为混合样本%s中没有其对应的reads。</font>\n' %(sample,files))
					skipped.append(sample)
					continue
				_=os.system('cat %s_f.2.fq %s_r.1.fq > %s_combined.2.fq' %(tmpbase+'fq/'+sample,tmpbase+'fq/'+sample,tmpbase+'fq/'+sample))
				_=os.system('rm %s_f.1.fq %s_r.2.fq %s_f.2.fq %s_r.1.fq' %(tmpbase+'fq/'+sample,tmpbase+'fq/'+sample,tmpbase+'fq/'+sample,tmpbase+'fq/'+sample))
				valid.append(j)
			data['file_name']=data['sample_id']+'_combined.1.fq;'+data['sample_id']+'_combined.2.fq'
			data=data.iloc[valid,]
			new_mixed=pd.concat([new_mixed,data])
		else:
			raise ValueError('\n## <font color=#FF0000>❌同一样本检测到2个以上文件，请检查。</font>\n')
	#生成新的mapping
	mapping=pd.concat([unmixed,new_mixed])

#若存在file相同、但没有index（非mixed）样本：
def find_duplicate_indices(lst):
    series = pd.Series(lst)
    duplicates = series[series.duplicated(keep=False)]
    duplicate_indices = duplicates.groupby(duplicates).apply(lambda x: list(x.index)).to_dict()
    return duplicate_indices

tmp = find_duplicate_indices(unmixed["file_name"])
pooled_ind = 0
if len(tmp) > 0:# 有pooled样本
	pooled_ind = 1
	pooled_index = [item for sublist in list(tmp.values()) for item in sublist]
	pooled = unmixed.iloc[pooled_index,]
	unmixed = unmixed.iloc[[i for i in list(range(unmixed.shape[0])) if i not in pooled_index],]
	mixed = mapping[mapping['mixed']=='yes']
	mapping = pd.concat([unmixed,mixed])

#对双端测序的进行拼接
def reverse_complement(seq):
	if len(seq)==0:
		return("")
	seq=str(seq)
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n', 'N': 'N'}
	bases = list(seq[::-1])
	letters = [complement[base] for base in bases]
	return(''.join(letters))

single=mapping[[(lambda x:True if x.count(';')==0 else False)(i) for i in mapping['file_name']]]
paired=mapping[[(lambda x:True if x.count(';')==1 else False)(i) for i in mapping['file_name']]]
if len(paired) > 0:
	print('\n# '+'开始对paired样本进行合并...'+'\n')
	new_name=[]
	for i in range(len(paired)):
		ind = 1
		sample=list(paired['sample_id'])[i]
		R1,R2=list(paired['file_name'])[i].split(';')
		R1=R1.strip('.gz')
		R2=R2.strip('.gz')
		
		print('\n## '+'开始对%s进行合并...' %(R1+';'+R2)+'\n')
		#如果已经合并，则跳过
		if os.path.exists(tmpbase+'fq/'+R1+'.merged.fq'):
			pass
			print('\n### '+'合并文件已存在，跳过合并...'+'\n')
		else:
			execut='vsearch --fastq_mergepairs %s --reverse %s --fastqout %s --fastq_allowmergestagger' %(tmpbase+'fq/'+R1,tmpbase+'fq/'+R2,tmpbase+'fq/'+R1+'.merged.fq')
			print(execut)
			val=os.system(execut)
			if val!=0:
				raise ValueError('\n### <font color=#FF0000>❌无法合并reads。</font>\n')
		
		#计算拼接后reads数比例
		before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+R1),shell=True).decode('utf-8').strip())
		after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+R1+'.merged.fq'),shell=True).decode('utf-8').strip())
		if before == 0:
			print('\n### <font color=#D89F12>⚠️警告：%s的reads数为0。</font>\n' %(R1+'、'+R2))
		elif after/before < 0.75:
			print('\n### <font color=#D89F12>⚠️警告：将%s合并后reads数从%d变为了%d，尝试使用单端分析...</font>\n' %(R1+'、'+R2,before/4,after/4))
			
			#若双端无法合并，则尝试使用单端reads
			amplicon=list(paired['amplicon_sequence'])[i].upper()
			sgRNA=list(paired['sgRNA_sequence'])[i].upper()
			#统计R1、R2的cover序列
			#R1
			end_length = 10
			R1_seqs=str(subprocess.check_output("seqkit sample --quiet -p %g %s|sed -n '2~4p'" %(min(100/(before/4),1), tmpbase+'fq/'+R1),shell=True).decode('utf-8').strip()).split('\n')
			R1_counter = Counter([ii[-end_length:] for ii in R1_seqs])
			R1_common_end = max(R1_counter, key=R1_counter.get)
			R2_seqs=str(subprocess.check_output("seqkit sample --quiet -p %g %s|sed -n '2~4p'" %(min(100/(before/4),1), tmpbase+'fq/'+R2),shell=True).decode('utf-8').strip()).split('\n')
			R2_counter = Counter([ii[-end_length:] for ii in R2_seqs])
			R2_common_end = max(R2_counter, key=R2_counter.get)
			if amplicon.find(R1_common_end) != -1 and reverse_complement(amplicon).find(R2_common_end) != -1:#R1为正向、R2为反向
				R1_seq = amplicon[:(amplicon.find(R1_common_end)+end_length)]
				R2_seq = reverse_complement(amplicon)[:(reverse_complement(amplicon).find(R2_common_end)+end_length)]
				if sgRNA in R1_seq:#选用R1
					print("\n### 选用R1进行分析(R1R2)\n")
					ind = 0
					new_name.append(R1)
					paired.iloc[i,paired.columns=='amplicon_sequence'] = R1_seq
					paired.iloc[i,paired.columns=='adapter_3'] = ''
					paired.iloc[i,paired.columns=='exclude_bp_from_right'] = 0
				elif reverse_complement(sgRNA) in R2_seq:#选用R2
					print("\n### 选用R2进行分析(R1R2)\n")
					ind = 0
					new_name.append(R2)
					paired.iloc[i,paired.columns=='amplicon_sequence'] = reverse_complement(R2_seq)
					paired.iloc[i,paired.columns=='adapter_5'] = list(paired['adapter_3'])[i]
					paired.iloc[i,paired.columns=='adapter_3'] = ''
					paired.iloc[i,paired.columns=='exclude_bp_from_left'] = 0
				else:
					print('\n### <font color=#D89F12>⚠️警告：无法采用单端分析，因为sgRNA不在R1或R2(R1R2)。</font>\n')
			elif reverse_complement(amplicon).find(R1_common_end) != -1 and amplicon.find(R2_common_end) != -1:#R1为反向、R2为正向
				R1_seq = reverse_complement(amplicon)[:(reverse_complement(amplicon).find(R1_common_end)+end_length)]
				R2_seq = amplicon[:(amplicon.find(R2_common_end)+end_length)]
				if reverse_complement(sgRNA) in R1_seq:#选用R1
					print("\n### 选用R1进行分析(R2R1)\n")
					ind = 0
					new_name.append(R1)
					paired.iloc[i,paired.columns=='amplicon_sequence'] = reverse_complement(R1_seq)
					paired.iloc[i,paired.columns=='adapter_3'] = ''
					paired.iloc[i,paired.columns=='exclude_bp_from_left'] = 0
				elif sgRNA in R2_seq:#选用R2
					print("\n### 选用R2进行分析(R2R1)\n")
					ind = 0
					new_name.append(R2)
					paired.iloc[i,paired.columns=='amplicon_sequence'] = R2_seq
					paired.iloc[i,paired.columns=='adapter_5'] = list(paired['adapter_3'])[i]
					paired.iloc[i,paired.columns=='adapter_3'] = ''
					paired.iloc[i,paired.columns=='exclude_bp_from_right'] = 0
				else:
					print('\n### <font color=#D89F12>⚠️警告：无法采用单端分析，因为sgRNA不在R1或R2(R2R1)。</font>\n')
			else:
				print("\n### 无法确定R1、R2方向，依旧使用双端分析。\n")
		else:
			print('完成。')
			#合并完成后删除R1和R2
			#_=os.system('rm %s %s' %(tmpbase+'fq/'+R1,tmpbase+'fq/'+R2))
		if ind == 1:
			new_name.append(R1+'.merged.fq')
		else:
			ind = 1
	paired['file_name']=new_name
mapping=pd.concat([single,paired])
if len(mapping) != 0:
	base=mapping[mapping['job']!='nnn']
	# base.to_csv(outbase+"mapping_base.csv",index=False)
	nnn=mapping[mapping['job']=='nnn']
else:
	base=mapping
	nnn=mapping
	
#分析base
test_control_lst=[]
test_control_list=[]
cut_mapping={}
if len(base)>0:
	print('\n# '+'开始对样本：%s进行base分析...' %(list(base['sample_id']))+'\n')
	for i in range(len(base)):
		sample=list(base['sample_id'])[i]
		file=list(base['file_name'])[i] 
		file=file.strip('.gz')
		stranded=list(base['stranded'])[i]
		sgRNA=list(base['sgRNA_sequence'])[i]
		amplicon=list(base['amplicon_sequence'])[i]
		adapter_5=list(base['adapter_5'])[i]
		adapter_3=list(base['adapter_3'])[i]
		mixed=list(base['mixed'])[i]
		paired=file.endswith('.merged.fq')
		quantification_center=list(base['quantification_center'])[i]
		if np.isnan(quantification_center):
			quantification_center=-3
		quantification_window_size=list(base['window_size'])[i]
		plot_window_size=list(base['window_size'])[i]
		if np.isnan(quantification_window_size):
			quantification_window_size=15
			plot_window_size=15
		exclude_bp_from_left=list(base['exclude_bp_from_left'])[i]
		if np.isnan(exclude_bp_from_left):
			exclude_bp_from_left=15
		exclude_bp_from_right=list(base['exclude_bp_from_right'])[i]
		if np.isnan(exclude_bp_from_right):
			exclude_bp_from_right=15
		alignment_cutoff=list(base['alignment_cutoff'])[i]
		if np.isnan(alignment_cutoff):
			alignment_cutoff=60
		control=list(base['control'])[i]
		if isinstance(control,str) and len(control)>0:
			if (sample,control) not in test_control_lst:
				test_control_lst.append((sample,control))
				test_control_list.append([sample,control,'',''])

		#检查是否为空，若是，则跳过：
		print('\n## '+'开始处理样本%s...' %(sample)+'\n')
		rows=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
		if rows==0:
			print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
			skipped.append(sample)
			continue
		if True:# mixed!='yes'
			if stranded=='yes':
				if isinstance(adapter_5,str) and len(adapter_5)>0:
					print('\n### '+'开始为%s去除5\'接头...' %(file)+'\n')
					val=os.system('cutadapt -g %s --discard-untrimmed -o %s %s -j 8' %(adapter_5,tmpbase+'fq/'+file+'.cut5.fq',tmpbase+'fq/'+file))
					if val!=0:
						raise ValueError('\n### <font color=#FF0000>❌为%s去除5\'接头时出错。</font>\n' %(file))
					#计算去接头后reads数比例
					before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
					after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut5.fq'),shell=True).decode('utf-8').strip())
					if after/before < 0.75:
						print('\n### <font color=#D89F12>⚠️警告：为%s去除5\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
					else:
						print('完成。')
					file=file+'.cut5.fq'
					if after==0:
						print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
						skipped.append(sample)
						continue
				if isinstance(adapter_3,str) and len(adapter_3)>0:
					#双端的3'接头需要进行反向互补
					if paired:
						adapter_3=reverse_complement(adapter_3)
					print('\n### '+'开始为%s去除3\'接头...' %(file)+'\n')
					val=os.system('cutadapt -a %s --discard-untrimmed -o %s %s -j 8' %(adapter_3,tmpbase+'fq/'+file+'.cut3.fq',tmpbase+'fq/'+file))
					if val!=0:
						raise ValueError('\n### <font color=#FF0000>❌为%s去除3\'接头时出错。</font>\n' %(file))
					#计算去接头后reads数比例
					before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
					after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut3.fq'),shell=True).decode('utf-8').strip())
					if after/before < 0.75:
						print('\n### <font color=#D89F12>⚠️警告：为%s去除3\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
					else:
						print('完成。')
					if isinstance(adapter_5,str) and len(adapter_5)>0:
						os.system('rm %s' %(tmpbase+'fq/'+file))
					file=file+'.cut3.fq'
					if after==0:
						print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
						skipped.append(sample)
						continue
			elif stranded!='yes':
				if isinstance(adapter_5,str) and len(adapter_5)>0:
					print('\n### '+'开始为%s去除5\'接头...' %(file)+'\n')
					val=os.system('cutadapt -g %s --revcomp --discard-untrimmed -o %s %s -j 8' %(adapter_5,tmpbase+'fq/'+file+'.cut5.fq',tmpbase+'fq/'+file))
					if val!=0:
						raise ValueError('\n### <font color=#FF0000>❌为%s去除5\'接头时出错。</font>\n' %(file))
					#计算去接头后reads数比例
					before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
					after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut5.fq'),shell=True).decode('utf-8').strip())
					if after/before < 0.75:
						print('\n### <font color=#D89F12>⚠️警告：为%s去除5\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
					else:
						print('完成。')
					file=file+'.cut5.fq'
					if after==0:
						print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
						skipped.append(sample)
						continue
					if isinstance(adapter_3,str) and len(adapter_3)>0:
						#双端的3'接头需要进行反向互补
						if paired:
							adapter_3=reverse_complement(adapter_3)
						print('\n### '+'开始为%s去除3\'接头...' %(file)+'\n')
						val=os.system('cutadapt -a %s --discard-untrimmed -o %s %s -j 8' %(adapter_3,tmpbase+'fq/'+file+'.cut3.fq',tmpbase+'fq/'+file))
						if val!=0:
							raise ValueError('\n### <font color=#FF0000>❌为%s去除3\'接头时出错。</font>\n' %(file))
						#计算去接头后reads数比例
						before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
						after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut3.fq'),shell=True).decode('utf-8').strip())
						if after/before < 0.75:
							print('\n### <font color=#D89F12>⚠️警告：为%s去除3\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
						else:
							print('完成。')
						#os.system('rm %s' %(tmpbase+'fq/'+file))
						file=file+'.cut3.fq'
						if after==0:
							print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
							skipped.append(sample)
							continue
				elif isinstance(adapter_3,str) and len(adapter_3)>0:
					#双端的3'接头需要进行反向互补
					if paired:
						adapter_3=reverse_complement(adapter_3)
					print('\n### '+'开始为%s去除3\'接头...' %(file)+'\n')
					val=os.system('cutadapt -a %s --revcomp --discard-untrimmed -o %s %s -j 8' %(adapter_3,tmpbase+'fq/'+file+'.cut3.fq',tmpbase+'fq/'+file))
					if val!=0:
						raise ValueError('\n### <font color=#FF0000>❌为%s去除3\'接头时出错。</font>\n' %(file))
					#计算去接头后reads数比例
					before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
					after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut3.fq'),shell=True).decode('utf-8').strip())
					if after/before < 0.75:
						print('\n### <font color=#D89F12>⚠️警告：为%s去除3\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
					else:
						print('完成。')
					file=file+'.cut3.fq'
					if after==0:
						print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
						skipped.append(sample)
						continue
		cut_mapping[sample]=file
		#调用CRISPResso2
		if isinstance(sgRNA,str) and len(sgRNA)>0:
			with_sgRNA=' -g '+sgRNA
		else:
			with_sgRNA=''
		# command='CRISPResso --max_rows_alleles_around_cut_to_plot 200 --name %s --amplicon_min_alignment_score %d --fastq_r1 %s --amplicon_seq %s%s -o %s --quantification_window_center %d --quantification_window_size %d --plot_window_size %d --exclude_bp_from_left %d --exclude_bp_from_right %d --allele_plot_pcts_only_for_assigned_reference --keep_intermediate --write_detailed_allele_table --min_frequency_alleles_around_cut_to_plot 0 -p 16' %(sample,alignment_cutoff,tmpbase+'fq/'+file,amplicon,with_sgRNA,outbase,quantification_center,quantification_window_size,plot_window_size,exclude_bp_from_left,exclude_bp_from_right)
		command='CRISPResso --max_rows_alleles_around_cut_to_plot 200 --name %s --amplicon_min_alignment_score %d --fastq_r1 %s --amplicon_seq %s%s -o %s --quantification_window_center %d --quantification_window_size %d --plot_window_size %d --exclude_bp_from_left %d --exclude_bp_from_right %d --allele_plot_pcts_only_for_assigned_reference --keep_intermediate --min_frequency_alleles_around_cut_to_plot 0 -p 16' %(sample,alignment_cutoff,tmpbase+'fq/'+file,amplicon,with_sgRNA,outbase,quantification_center,quantification_window_size,plot_window_size,exclude_bp_from_left,exclude_bp_from_right)
		print('\n### '+'开始进行处理...'+'\n')
		print(command)
		val=os.system(command)
		#tmp='CRISPResso_on_'+file.removesuffix('.fq').removesuffix('.fastq')
		#_=os.system('if [ -f '+tmp+'.html ];then mv '+tmp+'.html '+sample+'.html ;mv '+tmp+' '+sample+';fi')

def find_duplicate_indices(lst):
    index_dict = {}
    duplicate_indices = []
    for index, value in enumerate(lst):
        if value in index_dict:
            duplicate_indices.append(index)
        else:
            index_dict[value] = index
    return duplicate_indices


R1s = []
if pooled_ind:
	pooled_files=np.unique(pooled['file_name'])
	for i,files in enumerate(pooled_files):
		if files.count(';')==0:#为单端测序
			raise ValueError("单端测序尚未实现！")
		data=pooled[pooled['file_name']==files]
		R1,R2=files.split(';')
		R1=R1.strip('.gz')
		R2=R2.strip('.gz')
		f=open(tmpbase+'fq/'+R1+'.amplicon.txt','w')
		
		print('\n## '+'开始处理pooled样本（%s）：%s...' %(files, str(list(data['sample_id'])))+'\n')

		data_ind=0
		amplicon_set=set()
		dup_sample_set=set()
		for k in range(len(data)):
			sample=list(data['sample_id'])[k]
			sgRNA=list(data['sgRNA_sequence'])[k]
			amplicon=list(data['amplicon_sequence'])[k]
			control=list(data['control'])[k]
			if amplicon in amplicon_set:
				dup_sample_set.add(sample)
				continue
			else:
				amplicon_set.add(amplicon)
			if isinstance(control,str) and len(control)>0:
				if (sample,control) not in test_control_lst:
					test_control_lst.append((sample,control))
					test_control_list.append([sample,control,'CRISPRessoPooled_on_'+R1.replace('.','_')+'/',''])
			f.write(sample+"\t"+amplicon+"\t"+sgRNA+"\t"+"NA"+"\t"+"NA"+"\n")
			if data_ind==0:
				data_ind=1
				alignment_cutoff=list(data['alignment_cutoff'])[k]
				if np.isnan(alignment_cutoff):
					alignment_cutoff=60
				quantification_center=list(data['quantification_center'])[k]
				if np.isnan(quantification_center):
					quantification_center=-3
				quantification_window_size=list(data['window_size'])[k]
				plot_window_size=list(data['window_size'])[k]
				if np.isnan(quantification_window_size):
					quantification_window_size=15
					plot_window_size=15
				exclude_bp_from_left=list(data['exclude_bp_from_left'])[k]
				if np.isnan(exclude_bp_from_left):
					exclude_bp_from_left=15
				exclude_bp_from_right=list(data['exclude_bp_from_right'])[k]
				if np.isnan(exclude_bp_from_right):
					exclude_bp_from_right=15
		f.close()
		if len(dup_sample_set) > 0:
			print('\n### <font color=#D89F12>⚠️警告：'+f'跳过重复的样本:{list(dup_sample_set)}'+'</font>\n')
		command=f"CRISPRessoPooled -r1 {tmpbase+'fq/'+R1} -r2 {tmpbase+'fq/'+R2} -f {tmpbase+'fq/'+R1+'.amplicon.txt'} --max_rows_alleles_around_cut_to_plot 200 --name {R1.replace('.','_')} --amplicon_min_alignment_score {alignment_cutoff} -o {outbase} --quantification_window_center {quantification_center} --quantification_window_size {quantification_window_size} --plot_window_size {plot_window_size} --exclude_bp_from_left {exclude_bp_from_left} --exclude_bp_from_right {exclude_bp_from_right} --allele_plot_pcts_only_for_assigned_reference --keep_intermediate --write_detailed_allele_table --min_frequency_alleles_around_cut_to_plot 0 --min_reads_to_use_region 100 -p 16 --skip_failed"
		print('\n### '+'开始进行处理...'+'\n')
		print(command)
		val=os.system(command)
		if val == 0:
			R1s.append(R1.replace('.','_'))
	for i,files in enumerate(pooled_files):#重新遍历一次，更新control样本的路径prefix
		if files.count(';')==0:#为单端测序
			raise ValueError("单端测序尚未实现！")
		data=pooled[pooled['file_name']==files]
		R1,R2=files.split(';')
		R1=R1.strip('.gz')
		R2=R2.strip('.gz')
		amplicon_set=set()
		for k in range(len(data)):
			sample=list(data['sample_id'])[k]
			amplicon=list(data['amplicon_sequence'])[k]
			if amplicon in amplicon_set:
				continue
			else:
				amplicon_set.add(amplicon)
			if sample in [ss[1] for ss in test_control_lst]:
				control_ind = [index for index, value in enumerate(test_control_list) if value[1] == sample]
				for jk in control_ind:
					test_control_list[jk][3]='CRISPRessoPooled_on_'+R1.replace('.','_')+'/'

		
#summarize
if len(base)>0 or pooled_ind:
	print('\n# '+r'开始统计mapping statistics & modified% by indel...')
	_=os.system('/usr/bin/Rscript /var/www/platform/app/scripts/ngs/process_modified.R %s' %(" ".join([outbase]+[outbase+"CRISPRessoPooled_on_"+i for i in R1s])))

if 'test_control_list' in locals().keys():
	if len(test_control_list)>0:
		ff=open(outbase+"compare_summary.tsv",'w')
		ff.write("test\tcontrol\tSignificant_Deletions_site_num\tMaximum_efficiency_of_significant_Deletions\tMedian_efficiency_of_significant_Deletions\tSum_efficiency_of_significant_Deletions\tSignificant_Insertions_site_num\tMaximum_efficiency_of_significant_Insertions\tMedian_efficiency_of_significant_Insertions\tSum_efficiency_of_significant_Insertions\tSignificant_Substitutions_site_num\tMaximum_efficiency_of_significant_Substitutions\tMedian_efficiency_of_significant_Substitutions\tSum_efficiency_of_significant_Substitutions\tSignificant_Modifications_site_num\tMaximum_efficiency_of_significant_Modifications\tMedian_efficiency_of_significant_Modifications\tSum_efficiency_of_significant_Modifications\n")
		print('\n# '+'开始比较样本对：%s...' %(test_control_list)+'\n')
		for test,control,prefix_test,prefix_control in test_control_list:
			test = test.replace('.','_')
			control = control.replace('.','_')
			print('\n## '+'开始比较样本对%s...' %(str((test,control)))+'\n')
			ff.write(f"{test}\t{control}")
			if not os.path.exists(outbase+prefix_test+'CRISPResso_on_'+test):
				print('\n### <font color=#D89F12>⚠️警告：跳过，因为%s没有分析结果。</font>\n' %(test))
				ff.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
				continue
			if not os.path.exists(outbase+prefix_control+'CRISPResso_on_'+control):
				print('\n### <font color=#D89F12>⚠️警告：跳过，因为%s没有分析结果。</font>\n' %(control))
				ff.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
				continue
			command='CRISPRessoCompare -n %s -n1 %s -n2 %s -o %s %s %s' %(test+'_'+control,test,control,outbase,outbase+prefix_test+'CRISPResso_on_'+test,outbase+prefix_control+'CRISPResso_on_'+control)
			print(command)
			val=os.system(command)
			if val==0:
				wkdir=outbase+'CRISPRessoCompare_on_'+test+'_'+control+'/'
				o=open(wkdir+'Comparison_summary.txt','w')
				for cat in ['Deletions','Insertions','Substitutions','All_modifications']:
					file=cat+'_quantification.txt'
					values={}
					with open(wkdir+file,'r') as f:
						line=f.readline().strip().split()
						while line:
							values[line[0]]=line[1:]
							line = f.readline().strip().split()
					if sum([1 for i in values['qval_bonferroni'] if float(i)<0.05])==0:
						test_max=None
						test_median=None
						test_num=0
						test_sum=None
						ff.write("\t0\tNA\tNA\tNA")
					else:
						target=values[test+'_'+cat]
						target=[float(i) for i in target]
						total=values[test+'_total']
						total=[float(i) for i in total]
						q=values['qval_bonferroni']
						q=[float(i) for i in q]
						res=[target[i]/total[i] for i in range(len(target))]
						res=[abs(res[i]) for i in range(len(res)) if q[i]<0.05]
						test_num=len(res)
						test_max=max(res)
						test_median=np.median(res)
						test_sum=sum(res)

						o.write('#Significant '+cat+' sites:\n')
						o.write(str(test_num)+'\n\n')
						o.write('Maximum efficiency of significant ' + cat + ':\n')
						o.write(str(test_max)+'\n\n')
						o.write('Median efficiency of significant ' + cat + ':\n')
						o.write(str(test_median)+'\n\n')
						o.write('Sum(efficiency) of significant ' + cat + ':\n')
						o.write(str(test_sum)+'\n\n')

						ff.write(f"\t{test_num}\t{test_max}\t{test_median}\t{test_sum}")
				o.close()
				ff.write("\n")
			else:
				ff.write("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
		ff.close()

if len(nnn)>0:
	print('\n# '+'开始对样本：%s进行NNN分析...' %(list(nnn['sample_id']))+'\n')
	def gini(x):
		# (Warning: This is a concise implementation, but it is O(n**2)
		# in time and memory, where n = len(x).  *Don't* pass in huge
		# samples!)

		# Mean absolute difference
		mad = np.abs(np.subtract.outer(x, x)).mean()
		# Relative mean absolute difference
		rmad = mad/np.mean(x)
		# Gini coefficient
		g = 0.5 * rmad
		return g
	for i in range(len(nnn)):
		sample=list(nnn['sample_id'])[i]
		file=list(nnn['file_name'])[i]
		file=file.strip('.gz')
		stranded=list(nnn['stranded'])[i]
		amplicon=list(nnn['amplicon_sequence'])[i]
		adapter_5=list(nnn['adapter_5'])[i]
		adapter_3=list(nnn['adapter_3'])[i]
		mixed=list(nnn['mixed'])[i]
		paired=file.endswith('.merged.fq')
		print('\n## '+'开始处理样本%s...' %(sample)+'\n')

		#检查是否已经去过接头或者因为reads数为0被跳过
		if sample in skipped:
			continue
		if sample in cut_mapping.keys():
			file=cut_mapping[sample]
		else:
			#检查是否为空，若是，则跳过：
			rows=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
			if rows==0:
				print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
				skipped.append(sample+'-NNN')
				continue
			if stranded=='yes':
				if isinstance(adapter_5,str) and len(adapter_5)>0:
					print('\n### '+'开始为%s去除5\'接头...' %(file)+'\n')
					val=os.system('cutadapt -g %s --discard-untrimmed -o %s %s -j 8' %(adapter_5,tmpbase+'fq/'+file+'.cut5.fq',tmpbase+'fq/'+file))
					if val!=0:
						raise ValueError('\n### <font color=#FF0000>❌为%s去除5\'接头时出错。</font>\n' %(file))
					#计算去接头后reads数比例
					before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
					after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut5.fq'),shell=True).decode('utf-8').strip())
					if after/before < 0.75:
						print('\n### <font color=#D89F12>⚠️警告：为%s去除5\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
					else:
						print('完成。')
					file=file+'.cut5.fq'
					if after==0:
						print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
						skipped.append(sample+'-NNN')
						continue
				if isinstance(adapter_3,str) and len(adapter_3)>0:
					#双端的3'接头需要进行反向互补
					if paired:
						adapter_3=reverse_complement(adapter_3)
					print('\n### '+'开始为%s去除3\'接头...' %(file)+'\n')
					val=os.system('cutadapt -a %s --discard-untrimmed -o %s %s -j 8' %(adapter_3,tmpbase+'fq/'+file+'.cut3.fq',tmpbase+'fq/'+file))
					if val!=0:
						raise ValueError('\n### <font color=#FF0000>❌为%s去除3\'接头时出错。</font>\n' %(file))
					#计算去接头后reads数比例
					before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
					after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut3.fq'),shell=True).decode('utf-8').strip())
					if after/before < 0.75:
						print('\n### <font color=#D89F12>⚠️警告：为%s去除3\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
					else:
						print('完成。')
					if isinstance(adapter_5,str) and len(adapter_5)>0:
						os.system('rm %s' %(tmpbase+'fq/'+file))
					file=file+'.cut3.fq'
					if after==0:
						print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
						skipped.append(sample+'-NNN')
						continue
			elif stranded!='yes':
				if isinstance(adapter_5,str) and len(adapter_5)>0:
					print('\n### '+'开始为%s去除5\'接头...' %(file)+'\n')
					val=os.system('cutadapt -g %s --revcomp --discard-untrimmed -o %s %s -j 8' %(adapter_5,tmpbase+'fq/'+file+'.cut5.fq',tmpbase+'fq/'+file))
					if val!=0:
						raise ValueError('\n### <font color=#FF0000>❌为%s去除5\'接头时出错。</font>\n' %(file))
					#计算去接头后reads数比例
					before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
					after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut5.fq'),shell=True).decode('utf-8').strip())
					if after/before < 0.75:
						print('\n### <font color=#D89F12>⚠️警告：为%s去除5\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
					else:
						print('完成。')
					file=file+'.cut5.fq'
					if after==0:
						print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
						skipped.append(sample+'-NNN')
						continue
					if isinstance(adapter_3,str) and len(adapter_3)>0:
						#双端的3'接头需要进行反向互补
						if paired:
							adapter_3=reverse_complement(adapter_3)
						print('\n### '+'开始为%s去除3\'接头...' %(file)+'\n')
						val=os.system('cutadapt -a %s --discard-untrimmed -o %s %s -j 8' %(adapter_3,tmpbase+'fq/'+file+'.cut3.fq',tmpbase+'fq/'+file))
						if val!=0:
							raise ValueError('\n### <font color=#FF0000>❌为%s去除3\'接头时出错。</font>\n' %(file))
						#计算去接头后reads数比例
						before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
						after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut3.fq'),shell=True).decode('utf-8').strip())
						if after/before < 0.75:
							print('\n### <font color=#D89F12>⚠️警告：为%s去除3\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
						else:
							print('完成。')
						os.system('rm %s' %(tmpbase+'fq/'+file))
						file=file+'.cut3.fq'
						if after==0:
							print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
							skipped.append(sample+'-NNN')
							continue
				elif isinstance(adapter_3,str) and len(adapter_3)>0:
					#双端的3'接头需要进行反向互补
					if paired:
						adapter_3=reverse_complement(adapter_3)
					print('\n### '+'开始为%s去除3\'接头...' %(file)+'\n')
					val=os.system('cutadapt -a %s --revcomp --discard-untrimmed -o %s %s -j 8' %(adapter_3,tmpbase+'fq/'+file+'.cut3.fq',tmpbase+'fq/'+file))
					if val!=0:
						raise ValueError('\n### <font color=#FF0000>❌为%s去除3\'接头时出错。</font>\n' %(file))
					#计算去接头后reads数比例
					before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
					after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'.cut3.fq'),shell=True).decode('utf-8').strip())
					if after/before < 0.75:
						print('\n### <font color=#D89F12>⚠️警告：为%s去除3\'接头后reads数从%d变为了%d。</font>\n' %(file,before/4,after/4))
					else:
						print('完成。')
					file=file+'.cut3.fq'
					if after==0:
						print('\n### <font color=#D89F12>⚠️警告：跳过%s，因为对应的文件%s没有reads。</font>\n' %(sample,file))
						skipped.append(sample+'-NNN')
						continue
		#从amplicon_sequence提取NNN片段
		it = list(re.finditer("N+", amplicon.upper()))
		spans=[i.span() for i in it]
		print('NNN片段数：'+str(len(spans)))
		skip_merged_NNN=False
		if len(spans) > 1:
			dd=defaultdict(dict)
		for i,span in enumerate(spans):
			forw_seq = amplicon.upper()[:span[0]]
			back_seq = amplicon.upper()[span[1]:]
			#若forw_seq/back_seq中存在NNN，则相应延长该序列，使有效碱基数为17
			forw_seq = forw_seq[max(len(forw_seq)-forw_seq.count('N')-17,0):len(forw_seq)]
			back_seq = back_seq[0:(17+back_seq.count('N'))]
			print('\n### '+'开始为%s提取NNN序列%d...' %(file,i+1)+'\n')
			val=os.system('cutadapt -j 8 -g %s -o %s_inter_%d.fq %s -O 6 --discard-untrimmed --revcomp' %(forw_seq,tmpbase+'fq/'+file,i+1,tmpbase+'fq/'+file))
			if val!=0:
				raise ValueError('\n### <font color=#FF0000>❌为%s提取第%d个NNN序列时出错（5\'）。</font>\n' %(file,i+1))
				skip_merged_NNN=True
			#计算去5'后reads数比例
			before=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file),shell=True).decode('utf-8').strip())
			inter=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'_inter_'+str(i+1)+'.fq'),shell=True).decode('utf-8').strip())
			if inter/before < 0.75:
				print('\n### <font color=#D89F12>⚠️警告：为%s去除5\'序列后reads数从%d变为了%d。</font>\n' %(file,before/4,inter/4))
			else:
				print('完成。')
			if inter<4:
				print('\n### <font color=#D89F12>⚠️警告：跳过%s的第%s个NNN片段，因为对应的文件%s没有reads。</font>\n' %(sample,str(i+1),file+'_inter_'+str(i+1)+'.fq'))
				skipped.append(sample+'-NNN'+str(i+1))
				skip_merged_NNN=True
				_=os.system('rm %s' %(tmpbase+'fq/'+file+'_inter_'+str(i+1)+'.fq'))
				continue
			val=os.system('cutadapt -j 8 -a %s -o %s_count_%d.fq %s_inter_%d.fq -O 6 --discard-untrimmed;rm %s_inter_%d.fq' %(back_seq,tmpbase+'fq/'+file,i+1,tmpbase+'fq/'+file,i+1,tmpbase+'fq/'+file,i+1))
			if val!=0:
				raise ValueError('\n### 为%s提取第%d个NNN序列时出错（3\'）。\n' %(file,i+1))
			#计算去3'后reads数比例
			after=int(subprocess.check_output('wc -l %s|awk \'{print $1}\'' %(tmpbase+'fq/'+file+'_count_'+str(i+1)+'.fq'),shell=True).decode('utf-8').strip())
			if after/inter < 0.75:
				print('\n### <font color=#D89F12>⚠️警告：为%s去除3\'序列后reads数从%d变为了%d。</font>\n' %(file,inter/4,after/4))
			else:
				print('完成。')
			if after<4:
				print('\n### <font color=#D89F12>⚠️警告：跳过%s的第%s个NNN片段，因为对应的文件%s没有reads。</font>\n' %(sample,str(i+1),file+'_count_'+str(i+1)+'.fq'))
				skipped.append(sample+'-NNN'+str(i+1))
				skip_merged_NNN=True
				_=os.system('rm %s' %(tmpbase+'fq/'+file+'_count_'+str(i+1)+'.fq'))
				continue
			#NNN分析
			print('\n### 开始nnn分析...\n')
			input_file = tmpbase+'fq/'+file+'_count_'+str(i+1)+'.fq'
			seqs = SeqIO.parse(input_file, 'fastq')
			if len(spans) > 1:
				for record in seqs:
					dd[record.id]["s"+str(i+1)]=record.seq
				print("读取record数：",len(dd))
				seqs = SeqIO.parse(input_file, 'fastq')
			seq_freq = {}
			for seq in seqs:
				seq_freq.setdefault(str(seq.seq), 0)
				seq_freq[str(seq.seq)] += 1
			out_freq = tmpbase+'fq/'+sample+'_NNN_'+str(i+1)+'_freq.txt'
			# out_gini = outbase+sample+'_NNN_'+str(i+1)+'_gini.txt'
			# gini_index = gini(list(seq_freq.values()))
			# with open(out_gini, 'w') as g:
			#	 g.write(sample + '\t' + str(gini_index) + '\n')
			with open(out_freq, 'w') as f:
				for key, value in seq_freq.items():
					f.write(key + '\t' + str(value) + '\n')
			_=os.system('sort -k 2 -n -r %s > %s' %(out_freq,outbase+sample+'_NNN_'+str(i+1)+'_freq_sorted.txt'))
			# _=os.system('head -n 50 %s > %s' %(outbase+sample+'_NNN_'+str(i+1)+'_freq_sorted.txt',outbase+sample+'_NNN_'+str(i+1)+'_freq_sorted_top50.txt'))
			_=os.system('rm %s' %(out_freq))
			#_=os.system('rm %s' %(input_file))
			print('完成。')
		if len(spans) > 1 and skip_merged_NNN==False:
			print('\n### 开始merged nnn分析...\n')
			ddd=dict(dd)
			res={}
			com="+'|'+".join(["str(d.get('s"+str(i+1)+"'))" for i in range(len(spans))])
			for d in ddd.values():
				tmp = eval(com)
				res.setdefault(tmp, 0)
				res[tmp] += 1

			with open(tmpbase+'fq/'+sample+'_NNN_merged_freq.txt', 'w') as f:
				for key, value in res.items():
					f.write(key + '\t' + str(value) + '\n')
			_=os.system('sort -k 2 -n -r %s > %s' %(tmpbase+'fq/'+sample+'_NNN_merged_freq.txt',outbase+sample+'_NNN_merged_freq_sorted.txt'))
			_=os.system('rm %s' %(tmpbase+'fq/'+sample+'_NNN_merged_freq.txt'))
			print('完成。')
if len(skipped)>=1:
	print('\n# '+'注意：以下实例：%s没有分析结果，因为reads数为0。' %(str(skipped))+'\n')
print('\n# <font color=#00FF00>'+'✅所有分析完成。</font>'+'\n')

