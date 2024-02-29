(*This script computes the pass matrix from SH to HOT from L=0 
to a maximum desired maximum degree LMAX. It is meant to be
called from the command line like:
wolfram -rawterm < sh2hot.mathematica > sh2hothardcodes.h
 *)
LMAX = 16;
file = OpenWrite["sh2hothardcodes.cxx"];
WriteLine[file, "/**"];
WriteLine[file, 
  "* This file was auto-generated with Wolfram's Mathematica by \
running:"];
WriteLine[file, "*     $ wolfram -rawterm < sh2hot.mathematica"];
WriteString[file, "* In case you need SH maximum orders above L=", 
  LMAX, ", you can change the"];
WriteLine[file, ""];
WriteLine[file, 
  "* value of the variable LMAX in the very first line of"];
WriteLine[file, 
  "* sh2hot.mathematica and re-run it (it will take long!!!)"];
WriteLine[file, "*/"];
WriteLine[file, ""];
Print["#ifndef _sh2hothardcodes_cxx_"];
WriteLine[file, "#ifndef _sh2hothardcodes_cxx_"];
Print["#define _sh2hothardcodes_cxx_"];
WriteLine[file, "#define _sh2hothardcodes_cxx_"];
Print[""];
WriteLine[file, ""];
Print["#include \"sh2hot.h\""];
WriteLine[file, "#include \"sh2hot.h\""];
Print[""];
WriteLine[file, ""];
Print["namespace sh2hot"];
WriteLine[file, "namespace sh2hot"];
Print["{"];
WriteLine[file, "{"];
Print["    static const unsigned int SH2HOT_LMAX = ", LMAX, ";"];
WriteString[file, "    static const unsigned int SH2HOT_LMAX = ", 
  LMAX, ";"];
WriteLine[file, ""];
Print[""];
WriteLine[file, ""];
Print["    unsigned long sh2hotHardcodesDim( const unsigned int L \
)"];
WriteLine[file, 
  "    unsigned long sh2hotHardcodesDim( const unsigned int L )"];
Print["    {"];
WriteLine[file, "    {"];
Print["        unsigned long K = (L+1)*(L+2)/2;"];
WriteLine[file, "        unsigned long K = (L+1)*(L+2)/2;"];
Print["        return( L<=SH2HOT_LMAX ? K*K : 0 );"];
WriteLine[file, "        return( L<=SH2HOT_LMAX ? K*K : 0 );"];
Print["    }"];
WriteLine[file, "    }"];
Print[""];
WriteLine[file, ""];
Print["    void sh2hotHardcodes( const unsigned int L, BufferType buffer \
)"];
WriteLine[file, 
  "    void sh2hotHardcodes( const unsigned int L, BufferType buffer \
)"];
Print["    {"];
WriteLine[file, "    {"];
Print["        for( unsigned long k=0; k<sh2hotHardcodesDim(L); ++k \
)"];
WriteLine[file, 
  "        for( unsigned long k=0; k<sh2hotHardcodesDim(L); ++k )"];
Print["            buffer[k] = 0.0f;"];
WriteLine[file, "            buffer[k] = 0.0f;"];
Print[""];
WriteLine[file, ""];
Print["        switch(L){"];
WriteLine[file, "        switch(L){"];
For[L = 0, L <= LMAX, L += 2,
  K = (L + 1)*(L + 2)/2;
  Print["            case ", L, ":"];
  WriteString[file, "            case ", L, ":"];
  WriteLine[file, ""];
  r = 0;
  c = 0;
  For[l = 0, l <= L, l += 2,
   For[m = -l, m <= l, m++,
     If[m < 0,
      fa[z_] := (-1)^m*LegendreP[l, -m, z]*
        Sqrt[(2*l + 1)*Gamma[l + m + 1]/Gamma[l - m + 1]/(2*Pi)];
      fb[z_] := Cos[m*z];
      ,
      If[m > 0,
        fa[z_] := 
         LegendreP[l, m, z]*
          Sqrt[(2*l + 1)*Gamma[l - m + 1]/Gamma[l + m + 1]/(2*Pi)];
        fb[z_] := Sin[m*z];
        ,
        fa[z_] := LegendreP[l, 0, z]*Sqrt[(2*l + 1)/(4*Pi)];
        fb[z_] := 1;
        ];
      ];
     For[nx = L, nx >= 0, nx--,
      For[ny = L - nx, ny >= 0, ny--,
       nz = L - nx - ny;
       ia[w_] := fa[w]*w^nz*Sqrt[1 - w*w]^(nx + ny);
       ib[w_] := fb[w]*Cos[w]^nx*Sin[w]^ny;
       va = Integrate[ ia[w], {w, -1, 1} ];
       vb = Integrate[ ib[p], {p, 0, 2*Pi} ];
       res = 
        Simplify[ 
         va*vb*Gamma[L + 1]/Gamma[nx + 1]/Gamma[ny + 1]/
           Gamma[nz + 1] ];
       
       If[res == 0, NULL,
        Print["                buffer[", c, "*", K, "+", r, "] = ", 
         N[res, 19], ";"];
        WriteString[file, "                buffer[", c, "*", K, "+", 
         r, "] = ", N[res, 19], ";"];
        WriteLine[file, ""];
        ];
       c = c + 1;
       ]
      ];
     r = r + 1;
     c = 0;
     ];
   ];
  Print["                break;"];
  WriteLine[file, "                break;"];
  ];
Print["            default:"];
WriteLine[file, "            default:"];
Print["                break;"];
WriteLine[file, "                break;"];
Print["        } // End switch"];
WriteLine[file, "        } // End switch"];
Print["    } // End sh2hotHardcodes implementation"];
WriteLine[file, "    } // End sh2hotHardcodes implementation"];
Print[""];
WriteLine[file, ""];
Print["} // end namespace sh2hot"];
WriteLine[file, "} // end namespace sh2hot"];
Print[""];
WriteLine[file, ""];
Print["#endif"];
WriteLine[file, "#endif"];
Close[file];
