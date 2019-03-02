module AbelSand

output = include("Output.jl") # supports output functionality (e.g. array to csv)


const N_size = 700;    # the radius of the lattice Z^2, the actual size becomes (2*N+1)x(2*N+1)
const dx = [1,0,-1,0]; # for a given (x,y) in Z^2, (x + dx, y + dy) for all (dx,dy) covers the neighborhood of (x,y)
const dy = [0,1,0,-1];

struct L_coord
  # represents a lattice coordinate
  x::Int
  y::Int
end

function FindCoordinate(Z::Array{L_coord,1}, a::Int, b::Int)
   # in the given array Z of coordinates finds the (first) index of the tuple (a,b);
   # if no match, returns -1

  for i=1:length(Z)
    if (Z[i].x==a)&&(Z[i].y==b)
      return i;
    end
  end

  return -1;
end



function move(N)
    # the main function moving the pile sand grains of size N at the origin of Z^2 until the sandpile becomes stable

    Z_lat = zeros(UInt8, 2*N_size+1,2*N_size+1);  # models the integer lattice Z^2, we will have at most 4 sands on each vertex
    V_sites = falses(2*N_size+1,2*N_size+1);  # all sites which are visited by the sandpile process, are being marked here
    Odometer = zeros(UInt64, 2*N_size+1, 2*N_size+1); # stores the values of the odometer function


    walking = L_coord[];    # the coordinates of sites which need to move

    V_sites[N_size + 1,N_size + 1] = true;

    # i1, ... j2  -> show the boundaries of the box which is visited by the sandpile process
    i1, i2, j1, j2 = N_size + 1, N_size + 1, N_size + 1, N_size + 1 ;
    n = N

    t1 = time_ns();

    while n > 0
        n -= 1;

        Z_lat[N_size + 1, N_size + 1] += 1;
        if ( Z_lat[N_size + 1, N_size + 1] >= 4 )
            push!(walking, L_coord(N_size + 1, N_size + 1));
        end

        while ( length(walking) > 0 )
            w = pop!(walking);
            x = w.x;
            y = w.y;

            Z_lat[x, y] -= 4;
            Odometer[x,y] += 4;

            for k = 1:4
               Z_lat[x + dx[k], y + dy[k]] += 1; V_sites[x + dx[k], y + dy[k]] = true;
               if Z_lat[x + dx[k], y + dy[k]] >= 4
                  if FindCoordinate(walking, x + dx[k] , y + dy[k]) == -1
                      push!(walking, L_coord( x + dx[k], y + dy[k] ) );
                  end
               end
            end

            i1 = min(i1, x - 1);
            i2 = max(i2, x + 1);
            j1 = min(j1, y - 1);
            j2 = max(j2, y + 1);
        end

        #if n%200 == 0
        #    print(string("The boundaries are:: ", (i2-i1+1),"x",(j2-j1+1)), "\n");
        #    #flush(stdout);
        #end

    end #end of the main while
    t2 = time_ns();

    println("The final boundaries are:: ", (i2-i1+1),"x",(j2-j1+1), "\n");
    print( "time elapsed: " , (t2 - t1)/1.0e9, "\n" );

    Z_lat = output.printIntoFile(Z_lat,  0, string("Abel_Z_", N) );
    Odometer = output.printIntoFile(Odometer, 1,string("Abel_OD_", N)  );

    output.saveAsGrayImage(Z_lat, string("Abel_Z_", N), 20, 0);
    color_code = Dict(1=>[255, 128, 255], 2=>[255, 0, 0],3=>[0,128,255])
    output.saveAsRGBImage(Z_lat, string("Abel_Z_color_", N), color_code, 20, 0);

    # for the total elapsed time, it's better to use the @time macros on the main call


    return Z_lat, Odometer; # these are trimmed in output module

end # end of function move




end # end of the module
