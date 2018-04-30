#!/bin/bash

###########################################################################
#########							###########
#########		Autonenufaar				###########
######### @uthor : D Baux	david.baux<at>inserm.fr		###########
######### Date : 01/02/2017					###########
#########							###########
###########################################################################

###########################################################################
###########
########### 	Script to automate nenufaar pipeline
########### 	to treat NGS data
###########
###########################################################################


####	This script is meant to be croned
####	must check the runs directory, identify new runs
####	and launch nenufaar when a new run is available


##############		If any option is given, print help message	##################################
VERSION=1.0
USAGE="
Program: Autonenufaar
Version: ${VERSION}
Contact: Baux David <david.baux@inserm.fr>

Usage: This script is meant to be croned
	Should be executed once per 10 minutes

"


if [ $# -ne 0 ]; then
	echo "${USAGE}"
	echo "Error Message : Arguments provided"
	echo ""
	exit 1
fi


###############		Get options from conf file			##################################

CONFIG_FILE='/home/neuro_admin/autonenufaar/autonenufaar.conf'

#we check params against regexp

UNKNOWN=$(cat  ${CONFIG_FILE} | grep -Evi "^(#.*|[A-Z0-9_]*=[a-z0-9_ \.\/\$\{\}]*)$")
if [ -n "${UNKNOWN}" ]; then
	echo "Error in config file. Not allowed lines:"
	echo ${UNKNOWN}
	exit 1
fi

source ${CONFIG_FILE}

###############		1st check whether another instance of the script is running	##################

RESULT=$(ps x | grep -v grep | grep -c ${SERVICE})
#echo `ps x | grep -v grep |grep ${SERVICE} `
#echo "Result: ${RESULT}"

if [ "${RESULT}" -gt 3 ]; then
	exit 0
fi
#echo "Passed"

###############		Get run info file				 ##################################

# the file contains the run id and a code
# 0 => not treated => to do - used to retreat a run in case ex of error
# 1 => nenufaar is running -in case the security above does not work
# 2 => run treated - ignore directory
# the file is stored in an array and modified by the script

declare -A RUN_ARRAY #init array
while read LINE
do
	if echo ${LINE} | grep -E -v '^(#|$)' &>/dev/null; then
		if echo ${LINE} | grep -F '=' &>/dev/null; then
			RUN_ID=$(echo "${LINE}" | cut -d '=' -f 1)
			RUN_ARRAY[${RUN_ID}]=$(echo "${LINE}" | cut -d '=' -f 2-)
		fi
	fi
done < ${RUNS_FILE}


###############		Now we'll have a look at the content of the directories ###############################


#http://moinne.com/blog/ronald/bash/list-directory-names-in-bash-shell
#--time-style is used here to ensure awk $8 will return the right thing (dir name)
RUN_PATHS="${MINISEQ_RUNS_DIR} ${MISEQ_RUNS_DIR}"
for RUN_PATH in ${RUN_PATHS}
do
	RUNS=$(ls -l --time-style="long-iso" ${RUN_PATH} | egrep '^d' | awk '{print $8}' |  egrep '^[0-9]{6}_')
	for RUN in ${RUNS}
	do
		#echo ${RUN} ${RUN_ARRAY[${RUN}]}
		######do not look at runs set to 2 in the runs.txt file
		if [ -z "${RUN_ARRAY[${RUN}]}" ] || [ "${RUN_ARRAY[${RUN}]}" -eq 0 ]; then
		###### enable after testing
			#now we must look for the AnalysisLog.txt file
			#find "${RUN_PATH}${RUN}/" -mindepth 1 -maxdepth 3 -type f -name 'AnalysisLog.txt' -exec grep -e 'Total execution time' "{}" \;)
			#find ${RUN_PATH}${RUN} -mindepth 1 -maxdepth 3 -type f -name 'AnalysisLog.txt' -exec grep -e 'Total execution time' "{}" \;
			#get finished run
			if [ -n "$(find ${RUN_PATH}${RUN} -mindepth 1 -maxdepth 3 -type f -name 'AnalysisLog.txt' -exec grep -e 'Total execution time' '{}' \; -quit)" ]; then
				#get neurosensoriel run - manifest name must contain NS_
				#echo "${RUN_PATH}${RUN}/"
				for i in "${RUN_PATH}${RUN}/*.txt"; do
					#deal with multiple NS_ files
					j=0
					for FILE in ${i}; do
						if [[ "${FILE}" =~ /NS_.+_[[:digit:]]+\.txt ]] && [ "${j}" -eq 0 ]; then
							((j++))
							#we are in an NS run
							if [ -z "${RUN_ARRAY[${RUN}]}" ];then
								echo ${RUN}=1 >> ${RUNS_FILE}
								RUN_ARRAY[${RUN}]=1
							elif [ "${RUN_ARRAY[${RUN}]}" -eq 0 ];then
								#Change value on array and file to running
								sed -i -e "s/${RUN}=0/${RUN}=1/g" "${RUNS_FILE}"
								RUN_ARRAY[${RUN}]=1
							fi							
							#get fastqs in input
							#deal with manifest and intervals file - specific names choose yours
							INPUT='132_hg19'
							if [[ "${FILE}" =~ /NS_targeted_2\.txt ]]; then
								INPUT='2_hg19'
							elif [[ "${FILE}" =~ /NS_targeted_121_1\.txt ]]; then
								INPUT='121_hg19'
							fi
							#echo ${INPUT};exit;
							#LOG_FILE="${AUTONENUFAAR_DIR}autonenufaar.log"
							#touch ${LOG_FILE}
							#exec &>${LOG_FILE}
							echo "$(date) Copying fastqs of ${RUN} in input/NS/${INPUT}/${RUN} folder"
							mkdir "${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}"
							${RSYNC} -avq --exclude=L001 --exclude="*.txt" --exclude="Undetermined*" --exclude="*.xml" --exclude=Alignment --exclude=Matrix --exclude=Phasing --delete "${RUN_PATH}${RUN}/Data/Intensities/BaseCalls/" "${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}"
							#Then we need to organize it by reading the samplesheet
							#${GREP} -e '${SAMPLE_IDS}' "${RUNS_DIR}${RUN}/Alignment
							#f** it's a mess we'll do it directly with the fastqs
							mkdir "${RUN_PATH}${RUN}/nenufaar"
							chmod 777 "${RUN_PATH}${RUN}/nenufaar"
							
							#run CNV script
							#check if MobiCNV.py is present
							#https://github.com/mobidic/MobiCNV
							#requires xlsxwriter > 1.0.0 python module
							if [ -f "${MOBICNV}" ];then
								if [ -n "$(find ${MINISEQ_RUNS_DIR}${RUN}/Alignment_1/*/*_S1.coverage.csv -type f)" ];then
									echo "$(date) Running MobiCNV on run ${RUN}"
									echo "${PYTHON} ${MOBICNV} -i ${MINISEQ_RUNS_DIR}${RUN}/Alignment_*/*/ -t csv -o ${RUN_PATH}${RUN}/nenufaar/${RUN}.xlsx"
									${PYTHON} ${MOBICNV} -i ${MINISEQ_RUNS_DIR}${RUN}/Alignment_1/*/ -t csv -o ${RUN_PATH}${RUN}/nenufaar/${RUN}.xlsx
								elif [ -n "$(find ${MISEQ_RUNS_DIR}${RUN}/Data/Intensities/BaseCalls/*_S1.coverage.csv -type f)" ];then
									echo "$(date) Running MobiCNV on run ${RUN}"
									echo "${PYTHON} ${MOBICNV} -i ${MISEQ_RUNS_DIR}${RUN}/Data/Intensities/BaseCalls/ -t csv -o ${RUN_PATH}${RUN}/nenufaar/${RUN}.xlsx"
									${PYTHON} ${MOBICNV} -i ${MINISEQ_RUNS_DIR}${RUN}/Data/Intensities/BaseCalls/ -t csv -o ${RUN_PATH}${RUN}/nenufaar/${RUN}.xlsx
								fi
							fi
							
							SAMPLE_PATHS="$(ls ${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/*fastq.gz)"
							for SAMPLE_PATH in ${SAMPLE_PATHS}; do
								SAMPLE="$(basename "${SAMPLE_PATH}" | cut -d '_' -f 1)"
								if [ ! -d "${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/${SAMPLE}" ]; then
									mkdir "${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/${SAMPLE}"
									#echo "${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/${SAMPLE}_*fastq.gz ${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/${SAMPLE}/"
									mv ${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/${SAMPLE}_*fastq.gz ${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/${SAMPLE}/
									cp ${AUTONENUFAAR_DIR}empty_report.docx ${RUN_PATH}${RUN}/nenufaar/${SAMPLE}.docx
								fi
							done
							#$(ls input/NS/${INPUT}/${RUN}/*fastq.gz | cut -d '_' -f 1 | mkdir)

							#launch nenufaar
							echo "$(date) launching ${NENUFAAR} on run ${RUN}"
							cp ${NENUFAAR_DIR}input/NS/${INPUT}/*.list ${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/
							cp ${NENUFAAR_DIR}input/NS/${INPUT}/*.bed ${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/
							#${NOHUP} ${SH} ${NENUFAAR} -i "${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/" -o "${RUN_PATH}${RUN}/nenufaar/" -a annovar -f true -hsm true
							#mkdir -p "${NENUFAAR_DIR}output/NS/${RUN}"
							${NOHUP} ${SH} ${NENUFAAR} -i "${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/" -o "${NENUFAAR_DIR}output/NS/" -a annovar -f true -hsm true -l ns.txt -up false
							if [ "$?" -eq 0 ];then
								echo "$(date) Running MultiQC"
								${MULTIQC} -o "${NENUFAAR_DIR}output/NS/${RUN}" -n "${RUN}_multiqc" "${NENUFAAR_DIR}output/NS/${RUN}"
								echo "$(date) Moving ${RUN} in ${RUN_PATH}${RUN}/nenufaar/ folder"
								${RSYNC} -avq --remove-source-files "${NENUFAAR_DIR}output/NS/${RUN}" "${RUN_PATH}${RUN}/nenufaar/"
								#must be done form RS
								#chown -R LGM "${RUN_PATH}${RUN}/nenufaar/"
								rm -rf "${NENUFAAR_DIR}output/NS/${RUN}/"
								echo "$(date) Deleting local input files"
								rm -rf "${NENUFAAR_DIR}input/NS/${INPUT}/${RUN}/"
								##Change value on array and file to done
								sed -i -e "s/${RUN}=1/${RUN}=2/g" "${RUNS_FILE}"
								RUN_ARRAY[${RUN}]=2
								echo "$(date) ${NENUFAAR} done on run ${RUN}"

							else
								echo "$(date) ERROR in nenufaar execution: check log file in ${RUN_PATH}${RUN}/nenufaar/"
							fi

							#echo ${FILE}
						fi
					done
				done
			fi
		###### enable after testing
		fi
		###### enable after testing
	done
done

exit 0
