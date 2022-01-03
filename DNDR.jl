function DNDR(ALPHA,BETA)
      A=[-0.0180   -0.0280   -0.0370   -0.0480   -0.0430   -0.0520   -0.0620
         -0.0520   -0.0510   -0.0410   -0.0450   -0.0440   -0.0340   -0.0340
         -0.0520   -0.0430   -0.0380   -0.0450   -0.0410   -0.0360   -0.0270
         -0.0520   -0.0460   -0.0400   -0.0450   -0.0410   -0.0360   -0.0280
         -0.0540   -0.0450   -0.0400   -0.0440   -0.0400   -0.0350   -0.0270
         -0.0490   -0.0490   -0.0380   -0.0450   -0.0380   -0.0280   -0.0270
         -0.0590   -0.0570   -0.0370   -0.0470   -0.0340   -0.0240   -0.0230
         -0.0510   -0.0520   -0.0300   -0.0480   -0.0350   -0.0230   -0.0230
         -0.0300   -0.0300   -0.0270   -0.0490   -0.0350   -0.0200   -0.0190
         -0.0370   -0.0330   -0.0240   -0.0450   -0.0290   -0.0160   -0.0090
         -0.0260   -0.0300   -0.0190   -0.0330   -0.0220   -0.0100   -0.0250
         -0.0130   -0.0080   -0.0130   -0.0160   -0.0090   -0.0140   -0.0100]
      #
      S= 0.2 * ALPHA
      K= floor(S)+3
      if(K <= 1)
            K= 2
      end
      if(K >=  12)
            K=  11
      end
      DA= S - K +3
      L = K + fix( 1.1*sign(DA) )
      #
      S= 0.1 * BETA
      M= floor(S)+4
      if(M <= 1)
            M= 2
      end
      if(M >=  7)
            M=  6
      end
      DB= S - M + 4
      N= M + fix( 1.1*sign(DB) )
      #
      K=convert(Int64, K)
      L=convert(Int64, L)
      M=convert(Int64, M)
      N=convert(Int64, N)
      T= A[K,M]
      U= A[K,N]
      V= T + abs(DA) * (A[L,M] - T)
      W= U + abs(DA) * (A[L,N] - U)
      return dndrv = V + (W-V)  * abs(DB)
end
