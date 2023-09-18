#!/bin/bash
# -------------------- Time tracking, signal trapping, and requeue functions  ------------------------ 
secs2timestr() {
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "%02d:%02d:%02d\n" $h $m $s
}

timestr2secs() {
    echo $1| sed 's/-/:/' | awk -F: '{print $4, $3, $2, $1}'|awk '{print $1+60*$2+3600*$3+86400*$4}'
}

parse_job(){
    #set default
    #read <sig_time> from job script
    if [[ -z $ckpt_overhead ]]; then
        jscript=/var/spool/slurmd/job$SLURM_JOB_ID/slurm_script
        if [[ -f $jscript ]]; then
            x=`grep "#SBATCH*.--signal=" $jscript|grep -v "#*.#SBATCH"| tail -1 |awk -F@ '{print $2}'`
            if [[ -n $x ]]; then ckpt_overhead=`timestr2secs $x`; fi
        fi
    fi
    #read <comment> from job script
    if [[ -z $max_timelimit ]]; then
        jscript=/var/spool/slurmd/job$SLURM_JOB_ID/slurm_script
        if [[ -f $jscript ]]; then
            x=`grep -o '#SBATCH*.--time=[^ ]*' $jscript | awk -F'=' '{print $2}'`
            if [[ -n $x ]]; then max_timelimit=`timestr2secs $x`; fi

            x=`grep -o '#SBATCH*.-t=[^ ]*' $jscript | awk -F'=' '{print $2}'`
            if [[ -n $x ]]; then max_timelimit=`timestr2secs $x`; fi
        fi
    fi
    
    if [[ -z $ckpt_overhead ]]; then let ckpt_overhead=60; fi
    if [[ -z $max_timelimit ]]; then let max_timelimit=172800; fi

    echo checkpoint overhead \$ckpt_overhead: `secs2timestr $ckpt_overhead`
    echo maximum time limit  \$max_timelimit: `secs2timestr $max_timelimit`

    TOTAL_TIME=$(squeue -h -j $SLURM_JOB_ID -o %k)
    timeAlloc=$(squeue -h -j $SLURM_JOB_ID -o %l)

    fields=`echo $timeAlloc | awk -F ':' '{print NF}'`
    if [ $fields -le 2 ]; then
       timeAlloc=`echo 0:$timeAlloc`
    fi

    timeAlloc=`timestr2secs $timeAlloc`
    TOTAL_TIME=`timestr2secs $TOTAL_TIME`

    let remainingTimeSec=TOTAL_TIME-timeAlloc+ckpt_overhead
    if [ $remainingTimeSec -gt 0 ]; then
        remainingTime=`secs2timestr $remainingTimeSec`
        scontrol update JobId=$SLURM_JOB_ID Comment=$remainingTime

        let maxtime=`timestr2secs $max_timelimit`
        if [ $remainingTimeSec -gt $maxtime ]; then 
           requestTime=$max_timelimit
        else
           requestTime=$remainingTimeSec
        fi
        echo time remaining \$remainingTime: $remainingTime
        echo next timelimit \$requestTime: $requestTime
    fi
    requestTime=$((requestTime/60))        #convert to minutes instead of seconds
}

requeue_job() {

    parse_job 

    if [ -n $remainingTimeSec ] && [ $remainingTimeSec -gt 0 ]; then
        func="$1" ; shift
        for sig ; do
            trap "$func $sig" "$sig"
        done
    else 
       echo no more job requeues,done!
    fi
}

func_trap() {
######################################################
# -------------- checkpoint application --------------
######################################################
    # insert checkpoint command here if any
    $ckpt_command >&2
    #set -x
    trap '' SIGTERM
    scontrol requeue ${SLURM_JOB_ID} >&2
    scontrol update JobId=${SLURM_JOB_ID} TimeLimit=${requestTime} >&2
    trap - SIGTERM
    echo \$?: $? >&2
    #set +x
}

# Create dmtcp_command wrapper for easy communication with coordinator
dmtcp_command_job () {
    fname=dmtcp_command.$SLURM_JOB_ID	
    h=$1
    p=$2
    str="#!/bin/bash
    export PATH=$PATH
    export DMTCP_COORD_HOST=$h
    export DMTCP_COORD_PORT=$p
    dmtcp_command \$@"
    echo "$str" >$fname
    chmod a+rx $fname
}

restart_count () {
	echo ${SLURM_RESTART_COUNT:-0}
}


#----------------------------- Set up DMTCP environment for a job ------------#
start_coordinator()
{
    fname=dmtcp_command.$SLURM_JOBID
    h=`hostname`

    check_coordinator=`which dmtcp_coordinator`
    if [ -z "$check_coordinator" ]; then
        echo "No dmtcp_coordinator found. Check your DMTCP installation and PATH settings."
        exit 0
    fi

    #ZZdmtcp_coordinator --daemon --exit-on-last -p 0 --port-file $fname $@ 1>/dev/null 2>&1
    `which dmtcp_coordinator` --daemon --exit-on-last -p 0 --port-file $fname $@ 1>/dev/null 2>&1

    while true; do
        if [ -f "$fname" ]; then
            p=`cat $fname`
            if [ -n "$p" ]; then
                break
            fi
        fi
    done
    export DMTCP_COORD_HOST=$h
    export DMTCP_COORD_PORT=$p

    # Create dmtcp_command wrapper for easy communication with coordinator
    p=`cat $fname`
    str="#!/bin/bash
    export PATH=$PATH 
    export DMTCP_COORD_HOST=$h
    export DMTCP_COORD_PORT=$p
    dmtcp_command \$@"
    echo "$str" >$fname
    chmod a+rx $fname

    #log dmtcp
    echo $SLURM_JOB_ID $USER `date +%F` $(restart_count)  >> /usr/common/software/spool/dmtcp_command.log 
}

#wait for a process to complete
wait_pid () {
    pid=$1
    while [ -e /proc/$pid ]
    do
        sleep 1
    done
}

##function cr_run () {

###max_timelimit=00:02:30
###ckpt_overhead=30
###start_coordinator -i 30
##start_coordinator --exit-after-ckpt

##let restart_count=`squeue -h -O restartcnt -j $SLURM_JOB_ID`
##if (( restart_count == 0 )); then
##    dmtcp_launch --ckpt-signal 10 $@
##elif (( $restart_count > 0 )) && [[ -e dmtcp_restart_script.sh ]]; then
##    bash ./dmtcp_restart_script.sh -h $DMTCP_COORD_HOST -p $DMTCP_COORD_PORT
##else
##    echo "Failed to restart the job, exit"
##    exit
##fi
##}

#append_testpath () {
#    #export PATH=${PATH}:$DMTCP_DIR/test
#    export PATH=${PATH}:/global/common/sw/cray/cnl7/haswell/dmtcp/2019-10-24/test
#}

#prepend_testpath () {
#    #export PATH=$DMTCP_DIR/test:$PATH
#    export PATH=/global/common/sw/cray/cnl7/haswell/dmtcp/2019-10-24/test:$PATH
#}

#wait for coordinator to complete checkpointing
wait_coord () {
    let sum=0
    ckpt_done=0
    while true; do
        x=(`dmtcp_command.$SLURM_JOB_ID -s`)
        npeers=${x[6]/#*=/}
        running=${x[7]/#*=/}
        if [[ $npeers > 0 ]] && [[ $running == no ]] ; then
	    let sum=sum+1
            sleep 1
        elif [[ $npeers > 0 ]] && [[ $running == yes ]] ; then
	    ckpt_done=1
	    break
        else
	    break
	fi
    done
    if [[ $ckpt_done == 1 ]]; then
        echo checkpointing completed, overhead =  $sum seconds
    else
	echo no running job to checkpoint
    fi
}

#checkpoint before before requeue the job
ckpt_dmtcp () {
    dmtcp_command.$SLURM_JOB_ID -c 
    wait_coord
}

#added 9/6/2021 for MANA
wait_coord2 () {
    sleep 1 #in case there is a delay in changing the satus of the dmtcp_coordinator 
    let sum=0
    ckpt_done=0
    while true; do
        x=(`dmtcp_command -s`)
        #npeers=${x[6]/#*=/}
        #running=${x[7]/#*=/}
        npeers=${x[8]/#*=/}
        running=${x[9]/#*=/}
        if [[ $npeers > 0 ]] && [[ $running == no ]] ; then
            let sum=sum+1
            sleep 1
        elif [[ $npeers > 0 ]] && [[ $running == yes ]] ; then
            ckpt_done=1
            break
        else
            break
        fi
    done
    if [[ $ckpt_done == 1 ]]; then
        echo checkpointing completed, overhead =  $sum seconds
    else
        echo no running job to checkpoint
    fi
}


ckpt_mana () {
    dmtcp_command -c
    wait_coord2
}

