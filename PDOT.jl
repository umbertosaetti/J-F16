function PDOT(P3,P1)
    if (P1>=50.0)
        if (P3>=50.0)
          T=5.0
          P2=P1
        else
          P2=60.0
          T=RTAU(P2-P3)
        end
      else
        if (P3>=50.0)
          T=5.0
          P2=40.0
        else
          P2=P1
          T=RTAU(P2-P3)
        end
      end
      return pdotv=T*(P2-P3)
end
