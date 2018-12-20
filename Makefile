all: INPUT/vgrid.nc INPUT/topog.nc INPUT/forcing.nc

INPUT/vgrid.nc: indata.nc
	# vertical grid file for remapping
	ncks -O -v dz,zw,zt $< $@
	ncatted -a "units,z,c,c,m" $@

	#ncap2 -s 'zw=z; zt=z' -O $(TEMPNC) $@
	#@rm -f $(TEMPNC)

INPUT/topog.nc: indata.nc
	# topography
	ncks -O -v topog $< $@
	ncrename -v topog,depth $@

INPUT/forcing.nc: indata.nc
	$(eval TEMPNC := $(shell mktemp temp.XXX.nc))
	# extract prescribed fields: zonal wind stress sst restoring and salin flux
	ncks -O -v taux,sst,salin $< $(TEMPNC)
	ncrename -v salin,evap $(TEMPNC)

	# add new zeroed out forcing fields
	#   LW, SW, (evap), latent, sensible, liq_precip, froz_precip, liq_runoff,
	#   froz_runoff, SSS and tauy
	ncap2 -s 'LW[$$y,$$x]=0.0' \
		  -s 'SW[$$y,$$x]=0.0' \
		  -s 'latent[$$y,$$x]=0.0' \
		  -s 'sensible[$$y,$$x]=0.0' \
		  -s 'liq_precip[$$y,$$x]=0.0' \
		  -s 'froz_precip[$$y,$$x]=0.0' \
		  -s 'liq_runoff[$$y,$$x]=0.0' \
		  -s 'froz_runoff[$$y,$$x]=0.0' \
		  -s 'SSS[$$y,$$x]=35.0' \
		  -s 'tauy[$$y,$$x]=0.0' \
		  -s 'evap[$$y,$$x]=0.0' \
		  -O $(TEMPNC) $@

	@rm -f $(TEMPNC)

indata.nc: gendata.py
	# generate the initial data file
	python gendata.py
