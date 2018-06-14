function r = sobol_generate(m,n,skip)
% I4_SOBOL_GENERATE generates a Sobol dataset.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 December 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of points to generate.
%
%    Input, integer SKIP, the number of initial points to skip.
%
%    Output, real R(M,N), the points.
%
for j = 1 : n
  seed = skip + j - 1;
  [ r(1:m,j), seed ]  = i4_sobol ( m, seed );
end

end

% I4_SOBOL generates a new quasirandom Sobol vector with each call.
%
%  Discussion:
%
%    The routine adapts the ideas of Antonov and Saleev.
%
%    Thanks to Francis Dalaudier for pointing out that the range of allowed
%    values of DIM_NUM should start at 1, not 2!  17 February 2009.
%
%    This function was modified to use PERSISTENT variables rather than
%    GLOBAL variables, 13 December 2009.
%
%
%  Parameters:
%
%    Input, integer DIM_NUM, the number of spatial dimensions.
%    DIM_NUM must satisfy 1 <= DIM_NUM <= 40.
%
%    Input/output, integer SEED, the "seed" for the sequence.
%    This is essentially the index in the sequence of the quasirandom
%    value to be generated.  On output, SEED has been set to the
%    appropriate next value, usually simply SEED+1.
%    If SEED is less than 0 on input, it is treated as though it were 0.
%    An input value of 0 requests the first (0-th) element of the sequence.
%
%    Output, real QUASI(DIM_NUM), the next quasirandom vector.
%

%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 March 2012
%
%  Author:
%
%    Original FORTRAN77 version by Bennett Fox.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Antonov, Saleev,
%    USSR Computational Mathematics and Mathematical Physics,
%    Volume 19, 1980, pages 252 - 256.
%
%    Paul Bratley, Bennett Fox,
%    Algorithm 659:
%    Implementing Sobol's Quasirandom Sequence Generator,
%    ACM Transactions on Mathematical Software,
%    Volume 14, Number 1, pages 88-100, 1988.
%
%    Bennett Fox,
%    Algorithm 647:
%    Implementation and Relative Efficiency of Quasirandom 
%    Sequence Generators,
%    ACM Transactions on Mathematical Software,
%    Volume 12, Number 4, pages 362-376, 1986.
%
%    Ilya Sobol,
%    USSR Computational Mathematics and Mathematical Physics,
%    Volume 16, pages 236-242, 1977.
%
%    Ilya Sobol, Levitan, 
%    The Production of Points Uniformly Distributed in a Multidimensional 
%    Cube (in Russian),
%    Preprint IPM Akad. Nauk SSSR, 
%    Number 40, Moscow 1976.

function [quasi, seed] = i4_sobol (DimNum, seed)

persistent atmost;
persistent dim_max;
persistent dim_num_save;
persistent initialized;
persistent lastq;
persistent log_max;
persistent maxcol;
persistent poly;
persistent recipd;
persistent seed_save;
persistent v;

if ( isempty ( initialized ) )
  initialized = 0;
  dim_num_save = -1;
end

if nargin==1
  seed = tau_sobol(DimNum);
end

if ( ~initialized | DimNum ~= dim_num_save )

  initialized = 1;

  dim_max = 40;
  dim_num_save = -1;
  log_max = 30;
  seed_save = -1;
%
%  Initialize (part of) V.
%
  v(1:dim_max,1:log_max) = zeros(dim_max,log_max);

  v(1:40,1) = [ ...
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]';

  v(3:40,2) = [ ...
          1, 3, 1, 3, 1, 3, 3, 1, ...
    3, 1, 3, 1, 3, 1, 1, 3, 1, 3, ...
    1, 3, 1, 3, 3, 1, 3, 1, 3, 1, ...
    3, 1, 1, 3, 1, 3, 1, 3, 1, 3 ]';

  v(4:40,3) = [ ...
             7, 5, 1, 3, 3, 7, 5, ...
    5, 7, 7, 1, 3, 3, 7, 5, 1, 1, ...
    5, 3, 3, 1, 7, 5, 1, 3, 3, 7, ...
    5, 1, 1, 5, 7, 7, 5, 1, 3, 3 ]';

  v(6:40,4) = [ ...
                   1, 7, 9,13,11, ...
    1, 3, 7, 9, 5,13,13,11, 3,15, ...
    5, 3,15, 7, 9,13, 9, 1,11, 7, ...
    5,15, 1,15,11, 5, 3, 1, 7, 9 ]';

  v(8:40,5) = [ ...
                         9, 3,27, ...
   15,29,21,23,19,11,25, 7,13,17, ...
    1,25,29, 3,31,11, 5,23,27,19, ...
   21, 5, 1,17,13, 7,15, 9,31, 9 ]';

  v(14:40,6) = [ ...
            37,33, 7, 5,11,39,63, ...
   27,17,15,23,29, 3,21,13,31,25, ...
    9,49,33,19,29,11,19,27,15,25 ]';

  v(20:40,7) = [ ...
                                       13, ...
   33,115, 41, 79, 17, 29,119, 75, 73,105, ...
    7, 59, 65, 21,  3,113, 61, 89, 45,107 ]';

  v(38:40,8) = [ ...
                              7, 23, 39 ]';
%
%  Set POLY.
%
  poly(1:40)= [ ...
      1,   3,   7,  11,  13,  19,  25,  37,  59,  47, ...
     61,  55,  41,  67,  97,  91, 109, 103, 115, 131, ...
    193, 137, 145, 143, 241, 157, 185, 167, 229, 171, ...
    213, 191, 253, 203, 211, 239, 247, 285, 369, 299 ];

  atmost = 2^log_max - 1;
%
%  Find the number of bits in ATMOST.
%
  maxcol = i4_bit_hi1 ( atmost );
%
%  Initialize row 1 of V.
%
  v(1,1:maxcol) = 1;

end
%
%  Things to do only if the dimension changed.
%
if ( DimNum ~= dim_num_save )
%
%  Check parameters.
%
  if ( DimNum < 1 | dim_max < DimNum )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'I4_SOBOL - Fatal error!\n' );
    fprintf ( 1, '  The spatial dimension DIM_NUM should satisfy:\n' );
    fprintf ( 1, '    1 <= DIM_NUM <= %d\n', dim_max );
    fprintf ( 1, '  But this input value is DIM_NUM = %d\n', DimNum );
    return
  end

  dim_num_save = DimNum;
%
%  Initialize the remaining rows of V.
%
  for i = 2 : DimNum
%
%  The bits of the integer POLY(I) gives the form of polynomial I.
%
%  Find the degree of polynomial I from binary encoding.
%
    j = poly(i);
    m = 0;

    while ( 1 )

      j = floor ( j / 2 );

      if ( j <= 0 )
        break;
      end

      m = m + 1;

    end
%
%  Expand this bit pattern to separate components of the logical array INCLUD.
%
    j = poly(i);
    for k = m : -1 : 1
      j2 = floor ( j / 2 );
      includ(k) = ( j ~= 2 * j2 );
      j = j2;
    end
%
%  Calculate the remaining elements of row I as explained
%  in Bratley and Fox, section 2.
%
    for j = m + 1 : maxcol 
      newv = v(i,j-m);
      l = 1;
      for k = 1 : m
        l = 2 * l;
        if ( includ(k) )
          newv = bitxor ( newv, l * v(i,j-k) );
        end
      end
      v(i,j) = newv;
    end
  end
%
%  Multiply columns of V by appropriate power of 2.
%
  l = 1;
  for j = maxcol-1 : -1 : 1
    l = 2 * l;
    v(1:DimNum,j) = v(1:DimNum,j) * l;
  end
%
%  RECIPD is 1/(common denominator of the elements in V).
%
  recipd = 1.0 / ( 2 * l );

  lastq(1:DimNum) = 0;

end

seed = floor ( seed );

if ( seed < 0 )
  seed = 0;
end

if ( seed == 0 )

  l = 1;
  lastq(1:DimNum) = 0;

elseif ( seed == seed_save + 1 )
%
%  Find the position of the right-hand zero in SEED.
%
  l = i4_bit_lo0 ( seed );

elseif ( seed <= seed_save )

  seed_save = 0;
  l = 1;
  lastq(1:DimNum) = 0;

  for seed_temp = seed_save : seed - 1
    l = i4_bit_lo0 ( seed_temp );
    for i = 1 : DimNum
      lastq(i) = bitxor ( lastq(i), v(i,l) );
    end
  end

  l = i4_bit_lo0 ( seed );

elseif ( seed_save + 1 < seed )

  for seed_temp = seed_save + 1 : seed - 1
    l = i4_bit_lo0 ( seed_temp );
    for i = 1 : DimNum
      lastq(i) = bitxor ( lastq(i), v(i,l) );
    end
  end

  l = i4_bit_lo0 ( seed );

end
%
%  Check that the user is not calling too many times!
%
if ( maxcol < l )
  fprintf ( 1, '\n' );
  fprintf ( 1, 'I4_SOBOL - Fatal error!\n' );
  fprintf ( 1, '  Too many calls!\n' );
  fprintf ( 1, '  MAXCOL = %d\n', maxcol );
  fprintf ( 1, '  L =      %d\n', l );
  return
end
%
%  Calculate the new components of QUASI.
%
for i = 1 : DimNum
  quasi(i) = lastq(i) * recipd;
  lastq(i) = bitxor ( lastq(i), v(i,l) );
end

seed_save = seed;
seed = seed + 1;

end

% -------------------------------------------------------------------------
% Helper functions

function tau = tau_sobol(dim_num)
% TAU_SOBOL defines favorable starting seeds for Sobol sequences.
%
%  Discussion:
%
%    For spatial dimensions 1 through 13, this routine returns
%    a "favorable" value TAU by which an appropriate starting point
%    in the Sobol sequence can be determined.
%
%    These starting points have the form N = 2**K, where
%    for integration problems, it is desirable that
%      TAU + DIM_NUM - 1 <= K
%    while for optimization problems, it is desirable that
%      TAU < K.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 December 2009
%
%  Author:
%
%    Original FORTRAN77 version by Bennett Fox.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    IA Antonov, VM Saleev,
%    USSR Computational Mathematics and Mathematical Physics,
%    Volume 19, 1980, pages 252 - 256.
%
%    Paul Bratley, Bennett Fox,
%    Algorithm 659:
%    Implementing Sobol's Quasirandom Sequence Generator,
%    ACM Transactions on Mathematical Software,
%    Volume 14, Number 1, pages 88-100, 1988.
%
%    Bennett Fox,
%    Algorithm 647:
%    Implementation and Relative Efficiency of Quasirandom
%    Sequence Generators,
%    ACM Transactions on Mathematical Software,
%    Volume 12, Number 4, pages 362-376, 1986.
%
%    Stephen Joe, Frances Kuo
%    Remark on Algorithm 659:
%    Implementing Sobol's Quasirandom Sequence Generator,
%    ACM Transactions on Mathematical Software,
%    Volume 29, Number 1, pages 49-57, March 2003.
%
%    Ilya Sobol,
%    USSR Computational Mathematics and Mathematical Physics,
%    Volume 16, pages 236-242, 1977.
%
%    Ilya Sobol, YL Levitan,
%    The Production of Points Uniformly Distributed in a Multidimensional
%    Cube (in Russian),
%    Preprint IPM Akad. Nauk SSSR,
%    Number 40, Moscow 1976.
%
%  Parameters:
%
%    Input, integer DIM_NUM, the spatial dimension.  Only values
%    of 1 through 13 will result in useful responses.
%
%    Output, integer TAU, the value TAU.
%
dim_max = 13;

tau_table = [  0,  0,  1,  3,  5, ...
               8, 11, 15, 19, 23, ...
              27, 31, 35 ];

if ( 1 <= dim_num & dim_num <= dim_max )
  tau = tau_table(dim_num);
else
  tau = - 1;
end

end

function bit = i4_bit_hi1(n)
% I4_BIT_HI1 returns the position of the high 1 bit base 2 in an integer.
%
%  Example:
%
%       N    Binary     BIT
%    ----    --------  ----
%       0           0     0
%       1           1     1
%       2          10     2
%       3          11     2 
%       4         100     3
%       5         101     3
%       6         110     3
%       7         111     3
%       8        1000     4
%       9        1001     4
%      10        1010     4
%      11        1011     4
%      12        1100     4
%      13        1101     4
%      14        1110     4
%      15        1111     4
%      16       10000     5
%      17       10001     5
%    1023  1111111111    10
%    1024 10000000000    11
%    1025 10000000001    11
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 February 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the integer to be measured.
%    N should be nonnegative.  If N is nonpositive, the value will always be 0.
%
%    Output, integer BIT, the number of bits base 2.
%
i = floor ( n );
bit = 0;

while ( 1 )

  if ( i <= 0 )
    break;
  end

  bit = bit + 1;
  i = floor ( i / 2 );

end

end

function bit = i4_bit_lo0(n)
% I4_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
%
%  Example:
%
%       N    Binary     BIT
%    ----    --------  ----
%       0           0     1
%       1           1     2
%       2          10     1
%       3          11     3 
%       4         100     1
%       5         101     2
%       6         110     1
%       7         111     4
%       8        1000     1
%       9        1001     2
%      10        1010     1
%      11        1011     3
%      12        1100     1
%      13        1101     2
%      14        1110     1
%      15        1111     5
%      16       10000     1
%      17       10001     2
%    1023  1111111111     1
%    1024 10000000000     1
%    1025 10000000001     1
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 February 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the integer to be measured.
%    N should be nonnegative.
%
%    Output, integer BIT, the position of the low 1 bit.
%
bit = 0;
i = floor ( n );

while ( 1 )

  bit = bit + 1;
  i2 = floor ( i / 2 );

  if ( i == 2 * i2 )
    break;
  end

  i = i2;

end

end

