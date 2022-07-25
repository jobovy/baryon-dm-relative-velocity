NEW=_new

.PHONY: z11mass.ps z11lumfunc.ps z11mvmh.ps z11lumfunc-comps.ps \
	testMergerTree.sav testMergerTree2.sav

z11mass.ps: 
	$(PYTHON) plot_massfunction.py \
	z11mass.sav -o z11mass.ps --zrei=11 \
	z11massavg.sav --addavg

z0mass.ps: 
	$(PYTHON) plot_massfunction.py \
	z0mass.sav -o z0mass.ps --zrei=0 \
	z0massavg.sav --addavg

z11lumfunc$(NEW).ps:
	$(PYTHON) plot_lumfunc.py \
	z11mass.sav z11lumfunc$(NEW).sav -o z11lumfunc_new.ps --zrei=11 --fs=0.01 \
	z11lumfuncavg$(NEW).sav --addavg

z11lumfunc-comps$(NEW).ps:
	$(PYTHON) plot_lumfunccomponents.py \
	z11mass.sav z11lumfunc$(NEW).sav -o z11lumfunc-comps$(NEW).ps --zrei=11 --fs=0.01

z0lumfunc$(NEW).ps: plot_z0lumfunc.py
	$(PYTHON) plot_z0lumfunc.py \
	-o $@ z0lumfunc_avg_vbc0$(NEW).sav z0lumfunc_avg_vbc1$(NEW).sav \
	z0lumfunc_avg_vbc2$(NEW).sav z0lumfunc_avg_vbc3$(NEW).sav

sdsslumfunc$(NEW).ps: plot_sdsslumfunc.py
	$(PYTHON) plot_sdsslumfunc.py \
	-o $@ z0sdsslumfunc_avg_vbc0$(NEW).sav z0sdsslumfunc_avg_vbc1$(NEW).sav \
	z0sdsslumfunc_avg_vbc2$(NEW).sav z0sdsslumfunc_avg_vbc3$(NEW).sav

#Non-paper figures
z0lumfunc_massonly.ps: plot_z0lumfunc.py
	$(PYTHON) plot_z0lumfunc.py \
	-o $@ z0lumfunc_avg_vbc0.sav z0lumfunc_avg_vbc1_massonly.sav \
	z0lumfunc_avg_vbc2_massonly.sav z0lumfunc_avg_vbc3_massonly.sav

z0lumfunc_sfonly.ps: plot_z0lumfunc.py
	$(PYTHON) plot_z0lumfunc.py \
	-o $@ z0lumfunc_avg_vbc0.sav z0lumfunc_avg_vbc1_sfonly.sav \
	z0lumfunc_avg_vbc2_sfonly.sav z0lumfunc_avg_vbc3_sfonly.sav

z11mvmh.ps:
	$(PYTHON) plot_mvmh.py \
	z11mass.sav z11lumfunc.sav -o z11mvmh.ps --zrei=11 --fs=0.03

z0lumfunc_vbc0$(NEW).ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc$(NEW).sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc0_*.sav --vbc=0 \
	-o $@ --saveavg=z0lumfunc_avg_vbc0$(NEW).sav

z0mass_vbc0.ps: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t mass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc0_*.sav --vbc=0 \
	-o $@ --saveavg=z0mass_avg_vbc0.sav

z0massvsz11mass_vbc0.ps: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t z0massvsz11mass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav \
	../mergerTrees/mergerTree_vbc0_*.sav --vbc=0 \
	-o $@

mvmass_vbc0.png: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t mvmass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav \
	../mergerTrees/mergerTree_vbc0_*.sav --vbc=0 \
	-o $@

z0sdsslumfunc_vbc0$(NEW).ps:
	$(PYTHON)  plot_mergertree.py -t sdsslumfunc \
	--lumfuncsav=z11lumfunc$(NEW).sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc0_*.sav --vbc=0 \
	-o $@ --sdssselect=sdssselect61.sav \
	--saveavg=z0sdsslumfunc_avg_vbc0$(NEW).sav --fs=0.01

z0lumfunc_vbc1$(NEW).ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc$(NEW).sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc1_*.sav --vbc=1 \
	-o $@ --saveavg=z0lumfunc_avg_vbc1$(NEW).sav

z0mass_vbc1.ps: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t mass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc1_*.sav --vbc=1 \
	-o $@ --saveavg=z0mass_avg_vbc1.sav

z0lumfunc_vbc1_massonly.ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc1_*.sav --vbc=0 \
	-o $@ --saveavg=z0lumfunc_avg_vbc1_massonly.sav

z0lumfunc_vbc1_sfonly.ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc0_*.sav --vbc=1 \
	-o $@ --saveavg=z0lumfunc_avg_vbc1_sfonly.sav

z0sdsslumfunc_vbc1$(NEW).ps:
	$(PYTHON)  plot_mergertree.py -t sdsslumfunc \
	--lumfuncsav=z11lumfunc$(NEW).sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc1_*.sav --vbc=1 \
	-o $@ --sdssselect=sdssselect61.sav \
	--saveavg=z0sdsslumfunc_avg_vbc1$(NEW).sav --fs=0.01

z0massvsz11mass_vbc1.ps: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t z0massvsz11mass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav \
	../mergerTrees/mergerTree_vbc1_*.sav --vbc=1 \
	-o $@

mvmass_vbc1.png: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t mvmass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav \
	../mergerTrees/mergerTree_vbc1_*.sav --vbc=1 \
	-o $@

z0lumfunc_vbc2$(NEW).ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc$(NEW).sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc2_*.sav --vbc=2 \
	-o $@ --saveavg=z0lumfunc_avg_vbc2$(NEW).sav

z0mass_vbc2.ps: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t mass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc2_*.sav --vbc=2 \
	-o $@ --saveavg=z0mass_avg_vbc2.sav

z0lumfunc_vbc2_massonly.ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc2_*.sav --vbc=0 \
	-o $@ --saveavg=z0lumfunc_avg_vbc2_massonly.sav

z0lumfunc_vbc2_sfonly.ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc0_*.sav --vbc=2 \
	-o $@ --saveavg=z0lumfunc_avg_vbc2_sfonly.sav

z0sdsslumfunc_vbc2$(NEW).ps:
	$(PYTHON)  plot_mergertree.py -t sdsslumfunc \
	--lumfuncsav=z11lumfunc$(NEW).sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc2_*.sav --vbc=2 \
	-o $@ --sdssselect=sdssselect61.sav \
	--saveavg=z0sdsslumfunc_avg_vbc2$(NEW).sav --fs=0.01

z0massvsz11mass_vbc2.ps: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t z0massvsz11mass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav \
	../mergerTrees/mergerTree_vbc2_*.sav --vbc=2 \
	-o $@

mvmass_vbc2.png: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t mvmass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav \
	../mergerTrees/mergerTree_vbc2_*.sav --vbc=2 \
	-o $@

z0lumfunc_vbc3$(NEW).ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc$(NEW).sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc3_*.sav --vbc=3 \
	-o $@ --saveavg=z0lumfunc_avg_vbc3$(NEW).sav

z0mass_vbc3.ps: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t mass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc3_*.sav --vbc=3 \
	-o $@ --saveavg=z0mass_avg_vbc3.sav

z0lumfunc_vbc3_massonly.ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc3_*.sav --vbc=0 \
	-o $@ --saveavg=z0lumfunc_avg_vbc3_massonly.sav

z0lumfunc_vbc3_sfonly.ps:
	$(PYTHON)  plot_mergertree.py -t lumfunc \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc0_*.sav --vbc=3 \
	-o $@ --saveavg=z0lumfunc_avg_vbc3_sfonly.sav

z0sdsslumfunc_vbc3$(NEW).ps:
	$(PYTHON)  plot_mergertree.py -t sdsslumfunc \
	--lumfuncsav=z11lumfunc$(NEW).sav --masssav=z11mass.sav --addavg \
	../mergerTrees/mergerTree_vbc3_*.sav --vbc=3 \
	-o $@ --sdssselect=sdssselect61.sav \
	--saveavg=z0sdsslumfunc_avg_vbc3$(NEW).sav --fs=0.01

z0massvsz11mass_vbc3.ps: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t z0massvsz11mass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav \
	../mergerTrees/mergerTree_vbc3_*.sav --vbc=3 \
	-o $@

mvmass_vbc3.png: plot_mergertree.py
	$(PYTHON)  plot_mergertree.py -t mvmass \
	--lumfuncsav=z11lumfunc.sav --masssav=z11mass.sav \
	../mergerTrees/mergerTree_vbc3_*.sav --vbc=3 \
	-o $@

testMergerTree.sav:
	$(PYTHON) run_mergertree.py --seed=1 \
	$@ z0mass.sav --zrei=11 --mres=1000000.0 #=10.**6.
testMergerTree2.sav:
	$(PYTHON) run_mergertree.py --seed=2 \
	$@ z0mass.sav --zrei=11 --mres=1000000.0 #=10.**6.

#vbc=0 merger trees
../mergerTrees/mergerTree_vbc0_%.sav:
	$(PYTHON) run_mergertree.py \
	--seed=$(subst .sav,,$(subst ../mergerTrees/mergerTree_vbc0_,,$@)) \
	--vbc=0 \
	$@ z0mass.sav --zrei=11 --mres=1000000.0 #=10.**6.

#vbc=1 merger trees
../mergerTrees/mergerTree_vbc1_%.sav:
	$(PYTHON) run_mergertree.py \
	--seed=$(subst .sav,,$(subst ../mergerTrees/mergerTree_vbc1_,,$@)) \
	--vbc=1 \
	$@ z0mass.sav --zrei=11 --mres=1000000.0 #=10.**6.

#vbc=2 merger trees
../mergerTrees/mergerTree_vbc2_%.sav:
	$(PYTHON) run_mergertree.py \
	--seed=$(subst .sav,,$(subst ../mergerTrees/mergerTree_vbc2_,,$@)) \
	--vbc=2 \
	$@ z0mass.sav --zrei=11 --mres=1000000.0 #=10.**6.

#vbc=3 merger trees
../mergerTrees/mergerTree_vbc3_%.sav:
	$(PYTHON) run_mergertree.py \
	--seed=$(subst .sav,,$(subst ../mergerTrees/mergerTree_vbc3_,,$@)) \
	--vbc=3 \
	$@ z0mass.sav --zrei=11 --mres=1000000.0 #=10.**6.
