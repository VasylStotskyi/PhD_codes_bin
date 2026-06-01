#!/bin/sh
SHELL=/bin/sh
PATH=/sbin:/bin:/usr/sbin:/usr/bin
MAILTO=root

/usr/bin/rsync -ruvv -e "ssh -v -o 'ProxyJump vstotsky@ela.cscs.ch'" vstotsky@daint.cscs.ch:/scratch/snx3000/vstotsky/saponite/* /home/disk1/CSCS/
/usr/bin/rsync -ruv --exclude='*.cube' --exclude='*.kp' --exclude='*.wfn' --exclude='*.bak-*' --exclude='*.restart' --exclude='*.Hessian' --exclude='pwscf.xml' --exclude='pwscf.pdos_atm#*' --exclude='el_plot'  --exclude='mulliken.dat-1.mulliken' --exclude='core.*' --exclude='pwscf.wfc*' --exclude='pwscf.mix*' --exclude='pwscf.restart_scf*'  --exclude='pwscf.update'  --exclude='pwscf.restart_k*' --exclude='fort.*'  --exclude='pwscf.save' merlin-l-001:/data/user/stotskyi_v/* /home/disk1/Merlin6/
#cp -ur /home/disk1/CSCS/ /srv/
#cp -ur /home/disk1/Merlin6/ /srv/
