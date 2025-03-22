# Unfolding example on DY mass distribution

## quick start
### first setup
```bash
# -- at KNU
ssh -Y <your account name>@cms.knu.ac.kr

bash # -- from csh to bash
# -- just to use ROOT; you can put below lines in your ~/.bashrc
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_100 x86_64-centos7-gcc10-opt # -- could take some time

# -- under an arbitrary your working directory...
git clone git@github.com:KyeongPil-Lee/DYUnfold.git
cd DYUnfold
```

### after the first setup
```bash
bash
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_100 x86_64-centos7-gcc10-opt
cd DYUnfold
```

### Run codes
```bash
root -l -b -q produce_respM.cc
root -l -b -q do_unfold.cc
```