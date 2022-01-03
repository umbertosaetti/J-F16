function CZ(ALPHA,BETA,EL)
      A= [0.770;   0.241;  -0.100;  -0.416;  -0.731; -1.053; -1.366;
          -1.646; -1.917; -2.120; -2.248; -2.229]
      #
      S= 0.2 * ALPHA
      K= floor(S)+3
      if(K <= 1)
          K= 2
      end
      if(K >=  12)
          K=  11
      end
      DA= S - K + 3
      L = K + fix( 1.1*sign(DA) )
      K=convert(Int64, K)
      L=convert(Int64, L)
      S= A[K] + abs(DA) * (A[L] - A[K])
      return czv = S*(1-(BETA/57.3)^2) - .19*(EL/25.0)
end
