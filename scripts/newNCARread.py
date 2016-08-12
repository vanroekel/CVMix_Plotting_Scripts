def turbRead(basedir,testcasedir,rdir,nz,nscl,turbselect):
#
# 
# Read in turbulent statistics from NCAR LES model
#
# basedir = the base directory where the LES data resides
# testcasedir = the directory (inside base directory data is drawn from
# rdir = the individual restart directory if needed
# nz = number of points in vertical
# nscl = number of scalar quantities

    import numpy as np
    import os
#basedir='/Users/LVanRoekel/Desktop/LES/'
#start_time=1
#end_time=20001
#hisint=20
#interval=500

    if type(turbselect) is not list:
        print "turbselections must be in the form of a list"
        return None

    if basedir[len(basedir)-1]=='/':
        tempbdir=basedir
    else:
        tempbdir=basedir+'/'
    if testcasedir[len(testcasedir)-1]=='/':
        tempcdir=testcasedir
    else:
        tempcdir=testcasedir+'/'
    if rdir[len(rdir)-1]=='/':
        temprdir=rdir
    else:
        temprdir=rdir+'/'

    fdir=tempbdir+tempcdir+temprdir

    from operator import mul

#   nfiles=(end_time-start_time)/interval
#    ntinfile=interval/histint
#    ntimes=ntinfile*nfiles
    
    turbvals={'wwsb':1,'engz':2,'engsbz':3,'englez':4,'uxym':5,'vxym':6,'wxym':7,'txym':8,'divz':9,'utle':10,'utsb':11, \
              'vtle':12,'vtsb':13,'wtle':14,'wtsb':15,'wt_tot':16,'z':17,'zz':18,'shrz':19,'buyz':20,'triz':21,'uwsb':22, \
              'vwsb':23,'uwle':24,'vwle':25,'uw_tot':26,'vw_tot':27,'wcube':28,'wfour':29,'tcube':30,'ups':31,'vps':32, \
              'wps':33,'tps':34,'t_rprod':35,'t_wq':36,'t_wp':37,'t_tau':38,'t_tran':39,'t_buoy':40,'t_diss':41, \
              't_sprod':42,'tpress':43,'upress':44,'vpress':45,'wpressx':46,'wpressy':47,'SGSWT':48,'SGSUW':49, \
              'SGSVW':50,'wwT':51,'wwU':52,'wwV':53,'tpcp':54,'spcp':55,'wtDT':56,'uwDT':57,'vwDT':58,'wTdot':59}
    sizearr=[nz,nz+2,nz+2,nz,nz+2,nz+2,nz+2,[nz+2,nscl],nz+2,[nz,nscl],[nz,nscl],[nz,nscl],[nz,nscl],[nz,nscl],\
                    [nz,nscl],[nz,nscl],nz+2,nz+2,nz,nz,nz,nz,nz,nz,nz,nz,nz,nz,nz,[nz,nscl],nz,nz,nz,[nz,nscl],nz,nz,\
                    nz,nz,nz,nz,nz,nz,[nz,nscl],nz,nz,nz,nz,[nz,nscl],nz,nz,[nz,nscl],nz,nz,[nz,nscl-1],[nz,nscl-1] \
                    ,[nz,nscl],nz,nz,[nz,nscl]]    
    timearr=[]
    turbs2=[]
    for i in range(len(turbselect)):
        if turbselect[i] in turbvals.keys():
            turbs2.append(turbselect[i])

    nturbvals=len(turbs2)
            
    if nturbvals >= 1:

        sizeT1=[nz,nz+2,nz+2,nz,nz+2,nz+2,nz+2,[nz+2,nscl],nz+2,[nz,nscl],[nz,nscl],[nz,nscl],[nz,nscl],[nz,nscl],\
                    [nz,nscl],[nz,nscl],nz+2,nz+2,nz,nz,nz,nz,nz,nz,nz,nz,nz,nz,nz,[nz,nscl],nz,nz,nz,[nz,nscl],nz,nz,\
                    nz,nz,nz,nz,nz,nz,[nz,nscl],nz,nz,nz,nz,[nz,nscl],nz,nz,[nz,nscl],nz,nz,[nz,nscl-1],[nz,nscl-1] \
                    ,[nz,nscl],nz,nz,[nz,nscl]]

        turbT1={'wwsb':1,'engz':2,'engsbz':3,'englez':4,'uxym':5,'vxym':6,'wxym':7,'txym':8,'divz':9,'utle':10,'utsb':11, \
                    'vtle':12,'vtsb':13,'wtle':14,'wtsb':15,'wt_tot':16,'z':17,'zz':18,'shrz':19,'buyz':20,'triz':21,'uwsb':22, \
                    'vwsb':23,'uwle':24,'vwle':25,'uw_tot':26,'vw_tot':27,'wcube':28,'wfour':29,'tcube':30,'ups':31,'vps':32, \
                    'wps':33,'tps':34,'t_rprod':35,'t_wq':36,'t_wp':37,'t_tau':38,'t_tran':39,'t_buoy':40,'t_diss':41, \
                    't_sprod':42,'tpress':43,'upress':44,'vpress':45,'wpressx':46,'wpressy':47,'SGSWT':48,'SGSUW':49, \
                    'SGSVW':50,'wwT':51,'wwU':52,'wwV':53,'tpcp':54,'spcp':55,'wtDT':56,'uwDT':57,'vwDT':58,'wTdot':59}

        sumval1=0
        for ii in range(len(sizeT1)):
            if type(sizeT1[ii])==int:
                sumval1+=sizeT1[ii]
            else:
                sumval1+=reduce(mul,sizeT1[ii],1)
            
        sizeT2=[nz,nz+2,nz+2,nz,nz+2,nz+2,nz+2,[nz+2,nscl],nz+2,[nz,nscl],[nz,nscl],[nz,nscl],[nz,nscl],[nz,nscl],\
                    [nz,nscl],[nz,nscl],nz+2,nz+2,nz,nz,nz,nz,nz,nz,nz,nz,nz,nz,nz,[nz,nscl],nz,nz,nz,[nz,nscl],nz,nz,\
                    nz,nz,nz,nz,nz,nz,[nz,nscl],nz,nz,nz,nz,[nz,nscl],nz,nz,[nz,nscl],nz,nz,[nz,nscl-1],[nz,nscl],nz,nz]

        turbT2={'wwsb':1,'engz':2,'engsbz':3,'englez':4,'uxym':5,'vxym':6,'wxym':7,'txym':8,'divz':9,'utle':10,'utsb':11, \
                    'vtle':12,'vtsb':13,'wtle':14,'wtsb':15,'wt_tot':16,'z':17,'zz':18,'shrz':19,'buyz':20,'triz':21,'uwsb':22, \
                    'vwsb':23,'uwle':24,'vwle':25,'uw_tot':26,'vw_tot':27,'wcube':28,'wfour':29,'tcube':30,'ups':31,'vps':32, \
                    'wps':33,'tps':34,'t_rprod':35,'t_wq':36,'t_wp':37,'t_tau':38,'t_tran':39,'t_buoy':40,'t_diss':41, \
                    't_sprod':42,'tpress':43,'upress':44,'vpress':45,'wpressx':46,'wpressy':47,'SGSWT':48,'SGSUW':49, \
                    'SGSVW':50,'wwT':51,'wwU':52,'wwV':53,'tpcp':54,'wtDT':55,'uwDT':56,'vwDT':57}

        sumval2=0
        for ii in range(len(sizeT2)):
            if type(sizeT2[ii])==int:
                sumval2+=sizeT2[ii]
            else:
                sumval2+=reduce(mul,sizeT2[ii],1)
        
        sizeT3=[nz,nz+2,nz+2,nz,nz+2,nz+2,nz+2,[nz+2,nscl],nz+2,[nz,nscl],[nz,nscl],[nz,nscl],[nz,nscl],[nz,nscl],\
                    [nz,nscl],[nz,nscl],nz+2,nz+2,nz,nz,nz,nz,nz,nz,nz,nz,nz,nz,nz,[nz,nscl],nz,nz,nz,[nz,nscl],nz,nz,\
                    nz,nz,nz,nz,nz,nz]

        turbT3={'wwsb':1,'engz':2,'engsbz':3,'englez':4,'uxym':5,'vxym':6,'wxym':7,'txym':8,'divz':9,'utle':10,'utsb':11, \
                    'vtle':12,'vtsb':13,'wtle':14,'wtsb':15,'wt_tot':16,'z':17,'zz':18,'shrz':19,'buyz':20,'triz':21,'uwsb':22, \
                    'vwsb':23,'uwle':24,'vwle':25,'uw_tot':26,'vw_tot':27,'wcube':28,'wfour':29,'tcube':30,'ups':31,'vps':32, \
                    'wps':33,'tps':34,'t_rprod':35,'t_wq':36,'t_wp':37,'t_tau':38,'t_tran':39,'t_buoy':40,'t_diss':41, \
                    't_sprod':42}
        sumval3=0
        for ii in range(len(sizeT3)):
            if type(sizeT3[ii])==int:
                sumval3+=sizeT3[ii]
            else:
                sumval3+=reduce(mul,sizeT3[ii],1)
        
        flist = []
        for (dirpath, dirnames, filenames) in os.walk(fdir):
            flist.extend(filenames)
            break
    
        flist=sorted(flist)
        nfiles=len(flist)
        ascii=[]
        binary=[]
        for i in range(nfiles):
            if flist[i][len(flist[i])-4:len(flist[i])]=='ieee':
                binary.append(flist[i])
            else:
                ascii.append(flist[i])
    
        if len(binary)!=len(ascii):
            raise ValueError('not enough input files in directory')
        else:
            nfiles=len(binary)
        
        ntinfile=[]
        lenturbs=[]
        sumvalarr=[]
        
        
        
        for i in range(nfiles):
            fsize=os.stat(fdir+binary[i]).st_size

            diff1=fsize/(16.0*sumval1)-round(fsize/(16.0*sumval1))
            diff2=fsize/(16.0*sumval2)-round(fsize/(16.0*sumval2))
            diff3=fsize/(16.0*sumval3)-round(fsize/(16.0*sumval3))
            if diff1 == 0:
#                sizearr=sizeT1;
#                turbvals=turbT1.copy();
                lenturbs.append(len(turbT1))
                sumvalarr.append(sumval1);
            elif diff2 == 0:
#                sizearr=sizeT2;
#                turbvals=turbT2.copy();
                lenturbs.append(len(turbT2));
                sumvalarr.append(sumval2);
            elif diff3 == 0:
#                sizearr=sizeT3;
#                turbvals=turbT3.copy();
                lenturbs.append(len(turbT3));
                sumvalarr.append(sumval3);
            else:
                raise ValueError("corrupted hist file, remove "+binary[i]+" and "+ascii[i])
            ntinfile.append(int(fsize/(16.0*sumvalarr[i])))
          
        ntimes=sum(ntinfile)
        Lcopy=[None]*(nturbvals+1)
        svals=turbvals.copy() 
        L=[None]*len(svals)
       
        
        for i in turbvals.keys():
            if type(sizearr[turbvals[i]-1]) == list:
                svals[i]=[ntimes]+sizearr[turbvals[i]-1]
            else:
                svals[i]=[ntimes]+[sizearr[turbvals[i]-1]]
            L[turbvals[i]-1]=np.zeros(svals[i])
        
        t=0
        for i in range(nfiles): 
            f=open(fdir+binary[i],'rb')
            f2=open(fdir+ascii[i],'r')


    # Read in time values from the ascii history files
    #
            lencol=4
            counter=t
            for line in f2:
                line=line.strip()
                columns=line.split()
                if lencol==4:
                    timearr.append(columns[0])
                    counter+=1
                lencol=len(columns)
    #
    #
    #
    
            for j in range(ntinfile[i]):
                countv=0
                if lenturbs[i]==57:
                    stopval=lenturbs[i]+1
                    for ij in range(stopval):
                        if ij !=55:
                            s_arr=L[ij][t,:].shape
                            nvals=reduce(mul,s_arr,1)
                            tempvals=np.fromfile(f,dtype=np.dtype('f4'),count=nvals)
                            countv+=nvals
                            if len(s_arr) > 1:
                                tempvs=tempvals.reshape(s_arr,order='F')
                                L[ij][t,:]=tempvs
                            else:
                                L[ij][t,:]=tempvals
                        else:
                            ff=np.fromfile(f,dtype=np.dtype('f4'),count=nz)
                else:
                    stopval=lenturbs[i]
                
                    for ij in range(stopval):
                        s_arr=L[ij][t,:].shape
                        nvals=reduce(mul,s_arr,1)
                        tempvals=np.fromfile(f,dtype=np.dtype('f4'),count=nvals)
                        countv+=nvals
                        if len(s_arr) > 1:
                            tempvs=tempvals.reshape(s_arr,order='F')
                            L[ij][t,:]=tempvs
                        else:
                            L[ij][t,:]=tempvals
    
                ff=np.fromfile(f,np.dtype('f4'),count=3*sumvalarr[i])
                t+=1    
            f.close()    

        for i in range(nturbvals):
            Lcopy[i]=L[turbvals[turbselect[i]]-1]
        Lcopy[nturbvals]=timearr
        return Lcopy
    else:
        print "Error chosen turb profile '"'%s'"' not in dataset choose one of %s" % (turbselect,turbvals.keys())
    

