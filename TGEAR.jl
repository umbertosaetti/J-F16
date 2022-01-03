function TGEAR(THTL)
      if (THTL<=0.77)
        return tgearv = 64.94*THTL
      else
        return tgearv = 217.38*THTL-117.38
      end
end
