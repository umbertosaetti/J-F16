function ADC(VT,ALT)
   R0 = 2.377E-3
   TFAC = 1.0 - 0.703E-5*ALT
   T = 519.0*TFAC
   if ALT > 35000.0
      T= 390.0
   end
   RHO = R0*(TFAC^4.14)
   dum = 1.4*1716.3*T
   MACH = VT/sqrt(dum)
   QBAR = 0.5*RHO*VT*VT
   PS = 1715.0*RHO*T
return MACH, QBAR
end
