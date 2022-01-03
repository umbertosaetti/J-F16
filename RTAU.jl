function RTAU(DP)
      if (DP<=25.0)
          return rtauv=1.0
      elseif (DP>=50.0)
          return rtauv=0.1
      else
          return rtauv=1.9-.036*DP
      end
end
