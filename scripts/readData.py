def readData(modelName, file1, renameDict, avInterval, testcase, file2='null'):

	import xarray as xray
	import pandas as pd
	import numpy as np
	from mpas_xarray import preprocess_mpas
	from netCDF4 import Dataset
	import datetime

	import calcFluxes
	reload(calcFluxes)

	if modelName == 'MOM':
		if file2 == 'null':
			print 'error MOM requires two files prog.nc and visc.nc'
			return None
		else:
			dar1 = xray.open_dataset(file1).isel(xh=0,yh=0,xq=0,yq=0);
			dar2 = xray.open_dataset(file2).isel(xh=0,yh=0,xq=0,yq=0);

			start = pd.to_datetime(str(dar2.Time[0].values))
			second = pd.to_datetime(str(dar2.Time[1].values))
			stop = pd.to_datetime(str(dar2.Time[len(dar2.Time.values)-1].values))

			fh = Dataset(file2, mode='r')
			times = fh.variables['Time'][:]

			dtimes = [datetime.datetime(2000,1,1,0,0,0) + datetime.timedelta(days=x) for x in times];
			dar1.coords['Time'] = pd.to_datetime(dtimes)
			dar2.coords['Time'] = pd.to_datetime(dtimes)

			dar = dar1.merge(dar2)

			dar.rename(renameDict, inplace = True);

			nt = dar.temp.shape[0]
			nz = dar.temp.shape[1]

			zmid = np.zeros((nt,nz))
			ztop = np.zeros((nt,nz+1))
			for i in range(nt):
				zmid[i,:] = -dar.zm[:].values
				ztop[i,:] = -dar.zt[:].values

			dar['zmid'] = xray.DataArray(zmid, coords=[dar.Time, dar.zm], dims=['Time','zm'])
			dar['ztop'] = xray.DataArray(ztop, coords=[dar.Time, dar.zt], dims=['Time','zt'])

	elif modelName == 'POP':
		dar = xray.open_mfdataset(file1)
		fh = Dataset(file1,mode='r')
		times = fh.variables['time'][:]
		dar.rename({'zt':'zmid','time':'Time'}, inplace = True)
		dtimes = [datetime.datetime(2000,1,1,0,0,0) + datetime.timedelta(seconds=x) for x in times]
		dar.coords['Time'] = pd.to_datetime(dtimes)
	
		dar.rename(renameDict, inplace = True);
		dar['zmid'] = xray.DataArray(-dar.zmid.values, coords=[dar.Time, dar.zm], dims=['Time','zm'])
		dar['ztop'] = xray.DataArray(-dar.ztop.values, coords=[dar.Time, dar.zt], dims=['Time','zt'])

		#compute N2 for POP  assume linear eos with correct coefficients
		density = 1000.0*(1.0 - 2.E-4*(dar.temp.values - 5.0) + 8.E-4*(dar.salt.values - 35.))
		nt = dar.temp.shape[0]
		nz = dar.temp.shape[1]

		N2 = np.zeros((nt,nz+1))
		for i in range(nz-1):
			dden = density[:,i] - density[:,i+1]
			dz = dar.zmid[:,i].values - dar.zmid[:,i+1].values 
			N2[:,i+1] = 9.8106/1000.0*dden/dz

		dar['N2'] = xray.DataArray(N2, coords=[dar.Time, dar.zt], dims=['Time','zt'])

	elif modelName == 'MPAS':
		dar = xray.open_mfdataset(file1, preprocess=lambda x: preprocess_mpas(x, yearoffset=2000)).isel(nCells=0)
		dar.rename(renameDict, inplace = True);

		nz = dar.dims['zm']
		arr = np.array((dar.ztop[:,nz-1].values - dar.layerThickness[:,nz-1].values))
		ztemp = np.vstack((dar.ztop.values.T, arr)).T
		dar['ztop'] = xray.DataArray(ztemp, coords=[dar.Time, dar.zt], dims=['Time','zt'])

	nt = dar.dims['Time']
        nz = dar.dims['zm']
	if testcase == 'ConvectNew':
		ustar = 0.0*np.ones(len(dar.Time.values))
		tempsfcflux = -75./4.2E6*np.ones(len(dar.Time.values))
		saltsfcflux = 0.0*np.ones(len(dar.Time.values))
		buoyflux = 9.8*(2E-4*tempsfcflux - 8E-4*saltsfcflux)*np.ones(len(dar.Time.values))
    	if testcase == 'ConvectSaltWind' or testcase == 'ConvectSalt':
		ustar = 0.01*np.ones(len(dar.Time.values))
		tempsfcflux = -75./4.2E6*np.ones(len(dar.Time.values))
		saltsfcflux = 1.585E-8*35*np.ones(len(dar.Time.values))
		buoyflux = 9.8*(2E-4*tempsfcflux - 8E-4*saltsfcflux)
        
	if testcase == 'noMLdiurnal':
		ustar = 0.0*np.ones(len(dar.Time.values))
		tempsfcflux = -75/4.2E6
		saltsfcflux = 1.585E-8*35
		bFlux1 = 9.8*(2E-4*tempsfcflux - 8E-4*saltsfcflux)
		maxD = -np.pi*(bFlux1 / (9.8*2E-4) * 4.2E6)


		swflux = [maxD * max(np.cos(2.*np.pi*((((x-dar.Time.values[0]) / np.timedelta64(1, 's'))/86400) - 0.5)), 0) for x in dar.Time.values];
		if 'deep' in file1.lower():
			tempsfcflux = -75/4.2E6*np.ones(len(swflux)) + np.array(swflux)/4.2E6
		else:
			tempsfcflux = -75./4.2E6*np.ones(len(swflux)) + np.array(swflux)/4.2E6*(1.0 - (0.67*np.exp(-1.0) + 0.33*np.exp(-1.0/17.)))

		saltsfcflux = 1.585E-8*35*np.ones(len(swflux))
		buoyflux = 9.8*(2E-4*tempsfcflux - 8E-4*saltsfcflux)
	if testcase == 'AllForcingML':

		ustar = 0.01*np.ones(len(dar.Time.values))
		tempsfcflux = -75/4.2E6
		saltsfcflux = 1.585E-8*35
		bFlux1 = 9.8*(2E-4*tempsfcflux - 8E-4*saltsfcflux)
		maxD = -np.pi*(bFlux1 / (9.8*2E-4) * 4.2E6)

		swflux = [maxD * max(np.cos(2.*np.pi*((((x-dar.Time.values[0]) / np.timedelta64(1, 's'))/86400) - 0.5)), 0) for x in dar.Time.values];

		if 'deep' in file1.lower():
			tempsfcflux = -75/4.2E6*np.ones(len(swflux)) + np.array(swflux)/4.2E6
		else:
			tempsfcflux = -75./4.2E6*np.ones(len(swflux)) + np.array(swflux)/4.2E6*(1.0 - (0.67*np.exp(-1.0) + 0.33*np.exp(-1.0/17.)))

		saltsfcflux = 1.585E-8*35*np.ones(len(swflux))
		buoyflux = 9.8*(2E-4*tempsfcflux - 8E-4*saltsfcflux)
	if testcase == 'Shear':
		tempsfcflux = 0.0*np.ones(len(dar.Time.values))
		saltsfcflux = 0.0*np.ones(len(dar.Time.values))
		buoyflux = 0.0*np.ones(len(dar.Time.values))
		ustar = 0.1*np.ones(len(dar.Time.values))
        
	bldMinWB = np.zeros(nt);

	dar['nlt_temp'] = xray.DataArray((dar.nlt.values.T*tempsfcflux).T, coords=[dar.Time, dar.zt], dims=['Time','zt'])
	dar['nlt_salt'] = xray.DataArray((dar.nlt.values.T*saltsfcflux).T, coords=[dar.Time, dar.zt], dims=['Time','zt'])
	uw = calcFluxes.getFluxes(nt, dar.U.values, dar.kv.values, dar.zmid.values)
	vw = calcFluxes.getFluxes(nt, dar.V.values, dar.kv.values, dar.zmid.values)
	wt = calcFluxes.getFluxes(nt, dar.temp.values, dar.kh.values, dar.zmid.values, -dar.nlt_temp.values)
	ws = calcFluxes.getFluxes(nt, dar.salt.values, dar.kh.values, dar.zmid.values, -dar.nlt_salt.values)
	wt[:,0] = tempsfcflux
	ws[:,0] = saltsfcflux
	uw[:,0] = ustar**2

	sref = dar.salt[0,0].values

	b = 9.8*(2E-4*(dar.temp - 4.75) - 8E-4*(dar.salt - sref))

	wb = 9.8*(2E-4*wt - 8E-4*ws)
	for i in range(nt):
		spot = wb[i].argmax()
		bldMinWB[i] = dar.ztop[i,spot].values
    
	dar['uw'] = xray.DataArray(uw, coords=[dar.Time, dar.zt], dims=['Time','zt'])
	dar['vw'] = xray.DataArray(vw, coords=[dar.Time, dar.zt], dims=['Time','zt'])
	dar['wt'] = xray.DataArray(wt, coords=[dar.Time, dar.zt], dims=['Time','zt'])
	dar['ws'] = xray.DataArray(ws, coords=[dar.Time, dar.zt], dims=['Time','zt'])
	dar['wb'] = xray.DataArray(wb, coords=[dar.Time, dar.zt], dims=['Time','zt'])
	dar['bldMinWB'] = xray.DataArray(bldMinWB, coords=[dar.Time], dims=['Time'])
	dar['b'] = xray.DataArray(b, coords=[dar.Time, dar.zm], dims=['Time','zm'])
	dar['buoyFlux'] = xray.DataArray(buoyflux, coords=[dar.Time], dims=['Time'])
    	dar['tempsfcflux'] = xray.DataArray(tempsfcflux, coords=[dar.Time], dims=['Time'])
    	dar['saltsfcflux'] = xray.DataArray(saltsfcflux, coords=[dar.Time], dims=['Time'])
	dar['ustar'] = xray.DataArray(ustar, coords=[dar.Time], dims=['Time'])

	dsAv = dar.resample(avInterval,'Time')
	return dar, dsAv


