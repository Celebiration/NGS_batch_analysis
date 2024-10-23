#!/bin/bash
file=$1
username=$2
jobname=$3
db_name="NewPlatform"
whoami
set -e

trap "python -u /var/www/platform/app/scripts/ngs/taskfailed.py $jobname $username" EXIT

Base_dir='/var/www/platform/public'
# Base_dir=$PWD
base_dir=''$Base_dir'/tmp/ngs/'
work_dir='/var/www/platform/app/scripts/ngs/'
source /home/yhh/.bashrc
source activate base
cd ''$base_dir'/data/'
echo -e "\n# Contents\n"
echo -e "\n[toc]\n"

#检查传入文件是否为zip文件
if [ "${file##*.}"x = "zip"x ];then
	#获取上传文件名
	dir=${file##*/}
	dir=${dir%.*}
	#以上传文件名创建文件夹，并在此路径下进行解压
	cd $dir
	#在/var/www/platform/public/tmp/result/下创建结果文件夹
	result=$base_dir/result/$dir
	mkdir $result
	# conda env list
	echo -e "\n# 解压文件\n"
	/home/yhh/Downloads/7zz x $file
	# unzip $file

	#进入到压缩文件
	a=`ls`
	if [ ! -f mapping.tsv ] && [ ! -f mapping.txt ] && [ ! -f mapping.csv ]
	then
		cd `find . -maxdepth 1 -type d |tail -1`
	fi
	#获取压缩文件夹当前路径下所有文件名
	files=`ls`
	#获取.txt文档的所有数据
	dos2unix mapping.*
	
	#解压文件
	tmpp=`find . -type f -name "*.gz"`
	if [ -n "$tmpp" ];then
	echo -e "\n***********正在解压*************\n"
	gunzip fq/*.gz
	sleep 1s
	echo -e "\n***********解压完成*************\n"
	fi
	echo -e "\n正在处理数据....\n" 
	#调用python 脚本 第一个参数是文件的识别码，第二个参数是当前的路径,第三个参数是当前用户
	python -u ''$work_dir'/process.py' $dir $PWD $username
else
	#获取上传文件名
	dir=${file##*/}
	dir=${dir%.*}
	#以上传文件名创建文件夹，并在此路径下进行解压
	#mkdir $dir
	cd $dir
	#在/var/www/platform/public/tmp/result/下创建结果文件夹
	result=$base_dir/result/$dir
	mkdir -p $result
	# whoami
	# conda env list
	# unzip $file

	#进入到压缩文件
	a=`ls`
	if [ ! -f mapping.tsv ] && [ ! -f mapping.txt ] && [ ! -f mapping.csv ]
	then
		cd `find . -maxdepth 1 -type d |tail -1`
	fi
	#获取压缩文件夹当前路径下所有文件名
	files=`ls`
	#获取.txt文档的所有数据
	dos2unix mapping.*
	
	#解压文件
	tmpp=`find . -type f -name "*.gz"`
	if [ -n "$tmpp" ];then
	echo -e "\n***********正在解压*************\n"
	gunzip fq/*.gz
	sleep 1s
	echo -e "\n***********解压完成*************\n"
	fi
	echo -e "\n正在处理数据....\n" 
	#调用python 脚本 第一个参数是文件的识别码，第二个参数是当前的路径,第三个参数是当前用户
	python -u ''$work_dir'/process.py' $dir $PWD $username	
fi

#在/var/www/platform/public/tmp/result/下生成结果文件result并压缩
#输出下载路径
#zip -r $result
sed -i '/error:/I s/^/\n### <font color=#D89F12>⚠️/g;/error:/I s/$/<\/font>/g' ${base_dir}/data/${dir}/log.txt
cp ${base_dir}/data/${dir}/log.txt ${base_dir}/data/${dir}/log.md
sed -i 's/^--$//g' ${base_dir}/data/${dir}/log.md
#pandoc --pdf-engine=xelatex ${base_dir}/data/${dir}/log.md -o ${base_dir}/data/${dir}/log.pdf --toc -V CJKmainfont="SimSun" -V mainfont="Consolas" -H /home/fengchr/configs/pandoc_options.sty --wrap=preserve --top-level-division=chapter

cp ${base_dir}/data/${dir}/mapping.* ${result}/
cp ${base_dir}/data/${dir}/log.* ${result}/
cd $result
cd ..
echo $dir
zip -r "$dir.zip" $dir
echo "$result"

# 更新数据库
python ''$work_dir'/updateMysql.py' $dir $username
# python ''$work_dir'/updateMysql.py' $dir $username $db_name

trap - EXIT