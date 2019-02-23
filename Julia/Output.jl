module Output

# supports output functionality for the results of the sandpile symulations
# outputs the final grid in CSV format, as well as an image file

using CSV
using DataFrames
using Images


function TrimZeros(A)
  # given an array A trims any zero rows/columns from its borders
  # returns a 4 tuple of integers, i1, i2, j1, j2, where the trimmed array corresponds to A[i1:i2, j1:j2]
  # A can be either numeric or a boolean array

    i1, j1 = 1, 1;
    i2, j2 = size(A);

    zz = typeof(A[1,1])(0);   # comparison of a value takes into account the type as well

    # i1 is the first row which has non zero element
    for i = 1:size(A,1)
      q = false;
      for k = 1:size(A,2)
        if A[i,k] != zz
          q = true;
          #print("found non-zero at i1, k  = ", i1, " ", k , " value = ", A[i1, k], " \n");
          i1 = i;
          break
        end
      end

      if q == true
        break
      end
    end

    # i2 is the first from below row with non zero element
    for i in size(A,1):-1:1
      q = false;
      for k = 1:size(A,2)
        if A[i,k] != zz
          q = true;
          i2 = i;
          break
        end
      end

      if q == true
        break
      end
    end

    # j1 is the first column with non zero element

    for j = 1:size(A,2)
      q = false;
      for k = 1:size(A,1)
        if A[k,j] != zz
          j1 = j;
          q = true;
          break
        end
      end

      if q == true
        break
      end
    end

    # j2 is the last column with non zero element

    for j in size(A,2):-1:1
      q=false;
      for k=1:size(A,1)
        if A[k,j] != zz
          j2 = j;
          q=true;
          break
        end
      end

      if q==true
        break
      end
    end

    return i1, i2, j1, j2;

end

function addLayerofZeros(A, extraLayer)
  # adds layer of zeros from all corners to the given array A

  if extraLayer <= 0
    return A;
  end

  N, M = size(A);


  Z = zeros( typeof(A[1,1]), N + 2*extraLayer, M + 2*extraLayer );
  Z[(extraLayer+1):(N + extraLayer ), (extraLayer+1):(M + extraLayer ) ] = A;

  return Z;

end




function printIntoFile(A, extraLayer, strFileName)
  # exports a 2d matrix A into a csv file
  # @extraLayer is an integers adding layer of 0-s sorrounding the output matrix

  i1, i2, j1, j2  = TrimZeros( A );
  A = A[i1:i2, j1:j2];

  A = addLayerofZeros(A, extraLayer);

  CSV.write(string(strFileName,".csv"), DataFrame(A), writeheader = false);

  return A;

end


function Array_magnifier(A, cell_mag, border_mag)
  # A is the main array; @cell_mag is the magnifying size of the cell,
  # @border_mag is the magnifying size of the border between lattice cells

  # creates a new array where each cell of the original array A appears magnified by size = cell_mag


  total_factor = cell_mag + border_mag;

  A1 = zeros( typeof(A[1,1]), total_factor*size(A,1), total_factor*size(A,2) );

  for i = 1:size(A,1)
       for j = 1:size(A,2)
            for u = ((i-1)*total_factor +1):(i*total_factor)
            for v = ((j-1)*total_factor + 1):(j*total_factor)
               if ( ( u - (i-1)*total_factor <= cell_mag  ) && ( v - (j-1)*total_factor <= cell_mag ) )
                  A1[u, v] = A[i,j] ;
               end
            end
            end
      end
  end

  return A1;

end


end # end of the module
