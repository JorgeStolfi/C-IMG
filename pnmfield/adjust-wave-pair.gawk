#! /usr/bin/gawk -f
# Last edited on 2003-04-22 23:06:21 by stolfi

BEGIN {
  abort = -1;
  if (niter == "") { niter = 200; }
  split("", x);
  split("", y);
  split("", z); # Given function values
  coord_max = 0;
  n = 0; # Number of data points
}

(abort >= 0) { exit abort; }

/^ *([#]|$)/ { next; }

/^ *[0-9]+ +[0-9]+ +[0-9]+ +[0-9]+ +[0-9]+ +[\/] +[0-9]+ *$/ {
  x[n] = $1;
  y[n] = $2;
  z[n,0] = $3/$7;
  z[n,1] = $4/$7;
  z[n,2] = $5/$7;
  xa = x[n]; if (xa < 0) { xa = -xa; }
  ya = y[n]; if (ya < 0) { ya = -ya; }
  if (xa > coord_max) { coord_max = xa; }
  if (ya > coord_max) { coord_max = ya; }
  n++;
  next;
}

/./ { data_error("bad format"); }

END {
  if (abort >= 0) { exit abort; }

  # Wave parameters:
  nwaves = 2;       # Number of waves.
  dim = 2*nwaves+1; # Number of basis functions.

  split("", xper); # xper[k] is the x-period of wave [k].
  split("", yper); # yper[k] is the y-period of wave [k].
  split("", per2); # per2[k] is  xper[k]^2 + yper[k]^2.  
  split("", L);    # L[j] is the mean level of the comp [j].
  split("", C);    # C[k,j] is the amplitude of cos(r) for comp [j] of wave[k].
  split("", S);    # S[k,j] is the amplitude of sin(r) for comp [j] of wave[k].

  adjust_waves();
  print_waves();
}

function adjust_waves(  iter)
{
  
  # The best solution so far:
  split("", xper_best);
  split("", yper_best);
  split("", L_best);
  split("", C_best);
  split("", S_best);
  e_best = 9999999999;  # To force an update
  
  guess_initial_periods();
  for (iter = 0; iter < niter; iter++)
    { refine_waves(); }
  restore_best()
}

function guess_initial_periods( )
{
  xper[0] = 500; yper[0] = 0;
  xper[1] = 0; yper[1] = 500;
  compute_per2();
  per_min = coord_max/3;   # Minimum modulus of the period vectors. */
  per_max = coord_max*5;   # Maximum modulus of the period vectors. */
  per_w = 2*per_max;       # Uncertainty in the period vectors. */
  compute_optimal_amplitudes();
  compute_energy();
  update_best();
}

funtion update_best(   k,j)
{
  for (k = 0; k < nwaves; k++)
    { xper_best[k] = xper[k];
      yper_best[k] = yper[k];
    }
  for (j = 0; j < 3; j++)
    { L_best[j] = L[j];
      for (k = 0; k < nwaves; k++)
        { C_best[k,j] = C[k,j];
          S_best[k,j] = S[k,j];
        }
    }
  e_best = e;
}

funtion restore_best(   k,j)
{
  for (k = 0; k < nwaves; k++)
    { xper[k] = xper_best[k];
      yper[k] = yper_best[k];
    }
  compute_per2();
  for (j = 0; j < 3; j++)
    { L[j] = L_best[j];
      for (k = 0; k < nwaves; k++)
        { C[k,j] = C_best[k,j];
          S[k,j] = S_best[k,j];
        }
    }
  e = e_best;
}

function refine_waves(   magfactor)
{
  adjust_periods();
  compute_optimal_amplitudes();
  compute_energy();
  magfactor = exp(log(2.0)/(2*nwaves));
  if (e < e_best) 
    { update_best();
      per_w = per_w/magfactor;
      if (per_w < 1.0) { per_w = 1.0; }
    }
  else
    { restore_best(); 
      per_w = per_w*magfactor;
      if (per_w > 2*per_max) { per_w = 2*per_max; }
    }
}

function adjust_periods(  k,p2,p2_max,p2_min)
{
  p2_max = per_max*per_max;
  p2_min = per_min*per_min;
  for (k = 0; k < nwaves; k++)
    { do
        { xper[k] = adjust_coord(xper[k], per_w);
          yper[k] = adjust_coord(yper[k], per_w);
          p2 = xper[k]*xper[k] + yper[k]*yper[k];
          printf "." > "/dev/stderr";
        }
      while ((xper[k] < 0) || (p2 > p2_max) || (p2 < p2_min));
    }
  compute_per2();
}

function adjust_coord(x,w,   delta,xa)
{ 
  xa = (x >= 0 ? x : -x)
  delta = exp((2*rand() - 1)*log(1.0 + w/x));
  return x*delta;
}

function compute_optimal_amplitudes()
{
  compute_rigidity_matrix();
  compute_right_hand_side();
  solve_system();
  unpack_solution();
}

function compute_rigidity_matrix(   r1,r2,i,zb1,zb2,sum)
{
  # Computes the system's matrix M[i,j],
  # where i,j are basis indices.
  
  split("", M);
  for (r1 = 0; r1 < dim; r1++)
    { for (r2 = r1; r2 < dim; r2++)
        { sum = 0.0;
          for (i = 0; i < n; i++);
            { zb1 = eval_basis(r1, x[i], y[i]);
              zb2 = eval_basis(r2, x[i], y[i]);
              sum += zb1*zb2;
            }
          M[r1,r2] = sum;
          if (r1 != r2) { M[r2,r1] = sum; }
        }
    }
}

function compute_right_hand_side(   r,i,zb,zdata)
{
  # Computes the right-hand side matrix B[i,j] 
  # where i = basis index, j = color channel.

  split("", B);
  for (r = 0; r < dim; r++)
    { for (j = 0; j < 3; j++)
        { sum = 0.0;
          for (i = 0; i < n; i++);
            { zb = eval_basis(r, x[i], y[i]);
              zdata = z[i,j];
              sum += zb*zdata;
            }
          B[r,j] = sum;
        }
    }
}
  
function eval_basis(r,xp,yp,    k,alpha)
{
  if (r == 0) { return 1; }
  k = int((r-1)/2);
  alpha = twopi*(xper[k]*xp + yper[k]*yp)/per2[k];
  if (r % 2 == 1)
    { return cos(alpha); }
  else
    { return sin(alpha); }
}

function solve_system(  h,i,j,k,cm,sm,mm,wi,wk,ti,tk)
{
  i = 0;
  for (h = 0; h < dim; h++)
    { # Try to set M[i,h] = 1, and M[k,h] = 0 for all k > i:
      for (k = i+1; k < dim; k++)
        { # Rotate rows i and k so as to clear out M[k,i]:
          cm = M[i,h]; sm = M[k,h]; mm = cm*cm+sm*sm;
          if (mm != 0.0) 
            { cm /= mm; sm /= mm;
              M[i,h] = 1.0;
              M[k,h] = 0.0;
              for (j = h+1; j < dim; j++)
                { wi = M[i,j]; wk = M[k,j];
                  ti =  cm*wi + sm*wk;
                  tk = -sm*wi + cm*wk;
                  M[i,j] = ti; M[k,j] = tk;
                }
              for (j = 0; j < 3; j++)
                { wi = B[i,j]; wk = B[k,j];
                  ti =  cm*wi + sm*wk;
                  tk = -sm*wi + cm*wk;
                  B[i,j] = ti; B[k,j] = tk;
                }
            }
        }
      # If M[i,h] != 0, set M[k,h] = 0 for all k < i:
      if (M[i,h] != 0)
        { assert(M[i,h] == 1);
          for (k = 0; k < i; k++)
            { # Combine row i into row k so as to clear M[k,i]:
              sm = M[k,h];
              M[k,h] = 0.0;
              for (j = h+1; j < dim; j++)
                { M[k,j] = M[k,j] - sm*M[i,j]; }
              for (j = 0; j < 3; j++)
                { B[k,j] = B[k,j] - sm*B[i,j]; }
            }
          i++;
        }
    }
}
              
function unpack_solution(  h,i,j,solh)
{
  # Extract the coefficients L, C[k], S[k] from the 
  # M and B matrices, assuming that M has been 
  # diagonalized as much as possible.
  for (j = 0; j < 3; j++)
    { i = 0;
      for (h = 0; h < dim; h++)
        { if (M[i,h] == 0)
            { solh = 0.0; }
          else
            { assert(M[i,h] == 1)
              solh = B[i,j];
              i++;
            }
          if (h == 0)
            { L = solh; }
          else 
            { k = int((h-1)/2);
              if (h % 2 == 1)
                { C[k] = solh; }
              else
                { S[k] = solh; }
            }
        }
    }
}
                  
 
function data_error(msg)
{ printf "%d: **%s\n", FNR, msg > "/dev/stderr"; 
  abort = 1; exit abort;
}

function arg_error(msg)
{ printf "**%s\n", msg > "/dev/stderr"; 
  abort = 1; exit abort;
}
