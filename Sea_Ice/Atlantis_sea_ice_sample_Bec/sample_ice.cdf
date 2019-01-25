netcdf SETAS_VMPAIce_14Jan2016 {
dimensions:
	t = UNLIMITED ; // (1 currently)
	b = 11 ;
variables:
	double Ice_Class1(t, b) ;
		Ice_Class1:_FillValue = 0.25 ;
	double Ice_Class2(t, b) ;
		Ice_Class2:_FillValue = 0.25 ;
	double Ice_Class3(t, b) ;
		Ice_Class3:_FillValue = 0.25 ;
	double Ice_Class4(t, b) ;
		Ice_Class4:_FillValue = 0.25 ;
	double t(t) ;
		t:units = "seconds since 1970-01-01 00:00:00 +10" ;
		t:dt = 43200. ;
	double total_depth(t, b) ;
		total_depth:_FillValue = 0. ;

// global attributes:
		:title = "trivial" ;
		:geometry = "VMPA_setas.bgm" ;
		:parameters = "" ;
		:history = "Wed Aug  6 17:09:40 2014: ncks -d t,0,0 SETAS_VMPAIce.nc SETAS_VMPAIceSmall.nc" ;
		:NCO = "4.0.8" ;
data:

 Ice_Class1 =
  _, 0.1, _, _, _, _, _, _, _, _, _ ;

 Ice_Class2 =
  _, 0.1, _, _, _, _, _, _, _, _, _ ;

 Ice_Class3 =
  _, 0.1, _, _, _, _, _, _, _, _, _ ;

 Ice_Class4 =
  _, 0.1, _, _, _, _, _, _, _, _, _ ;

 t = 43200 ;

 total_depth =
  _, 10, _, _, _, _, _, _, _, _, _ ;
}
