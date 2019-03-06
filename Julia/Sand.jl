# ==================================================
# Author: Hayk Aleksanyan
# ==================================================
module Sand

output = include("Output.jl") # supports output functionality (e.g. array to csv)


const N_size = 700;    # the radius of the lattice Z^2, the actual size becomes (2*N+1)x(2*N+1)
const dx = [1,0,-1,0]; # for a given (x,y) in Z^2, (x + dx, y + dy) for all (dx,dy) covers the neighborhood of (x,y)
const dy = [0,1,0,-1];

struct L_coord
  # represents a lattice coordinate
  x::Int
  y::Int
end

struct Graph_rep
  # used in representing a Graph, where each neigh array shows the indices of neighboring vertices to the given vertex
  neigh::Array{Int,1}
end

function Init_CheckArray()
  # returns an extended neighborhood of the origin of Z^2 lattice

  v=zeros(Int,(13,2));

  v[1,1]=0;     v[1, 2] = 0;
  v[2, 1] = 0;  v[2, 2] = 1;
  v[3, 1] = 0;  v[3, 2] = 2;
  v[4, 1] = 0;  v[4, 2] = -1;
  v[5, 1] = 0;  v[5, 2] = -2;
  v[6, 1] = 1;  v[6, 2] = 0;
  v[7, 1] = 2;  v[7, 2] = 0;
  v[8, 1] = -1; v[8, 2] = 0;
  v[9, 1] = -2; v[9, 2] = 0;
  v[10, 1] = 1;  v[10, 2] = 1;
  v[11, 1] = -1; v[11, 2] = 1;
  v[12, 1] = -1; v[12, 2] = -1;
  v[13, 1] = 1;  v[13, 2] = -1;

  return v;

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



function moveSand(T_mass)
    # the main function moving a sand of mass T_mass out of the origin of Z^2
    # implements Algorthim 2 of https://arxiv.org/abs/1607.01525 Discrete balayage and boundary sandpile
    # provides a precise compuation of the final distribution of the boundary sandpile

    Z_lat = zeros(2*N_size+1,2*N_size+1);     # models the integer lattice Z^2
    V_sites = falses(2*N_size+1,2*N_size+1);  # all sites which are visited by the sandpile process, are being marked here
    Odometer = zeros(2*N_size+1, 2*N_size+1); # stores the values of the odometer function


    w_Bound = L_coord[];  # the coordinates of boundary sites which need to move
    w_Int = L_coord[];    # the coordinates of interior sites which need to move

    v_check = Init_CheckArray();   #Array(Int,13,2);

    Norm_mass=T_mass^0.5;  # the capacity of the boundary

    Z_lat[N_size+1,N_size+1] = T_mass;  # the initial mass is stored at the origin
    V_sites[N_size+1,N_size+1] = true;

    push!(w_Bound,L_coord(N_size+1,N_size+1));  # the origin needs to topple

    i, j = 1, 0 ;

    # i1, ... j2  -> show the boundaries of the box which is visited by the sandpile process
    i1, i2, j1, j2 = N_size + 1, N_size + 1, N_size + 1, N_size + 1 ;

    iter_ = 0;   # number of total iterations

    while 0 == 0
        # we iterate until the process becomes stable
        # it is proved in our paper https://arxiv.org/abs/1607.01525
        # that the process stabilises after finitely many moves
        top_ = length(w_Bound);

        iter_ = iter_+1;
        println("\nIteration N:: ", iter_);
        println("The boundaries are:: ", (i2-i1+1),"x",(j2-j1+1), "\n");
        println("Moving the boundary...");

        t1 = time_ns();

        while (top_>0)
            x=w_Bound[top_].x;
            y=w_Bound[top_].y;

            # any mass at the (x,y) needs to be moved out
            dist_mass = Z_lat[x,y]/4;    # the mass is split evenly among lattice neighbors

            Z_lat[x,y + 1] = Z_lat[x,y + 1] + dist_mass; V_sites[x,y + 1] = true;
            Z_lat[x,y - 1] = Z_lat[x,y - 1] + dist_mass; V_sites[x,y - 1] = true;
            Z_lat[x - 1,y] = Z_lat[x - 1,y] + dist_mass; V_sites[x - 1,y] = true;
            Z_lat[x + 1,y] = Z_lat[x + 1,y] + dist_mass; V_sites[x + 1,y] = true;

            Odometer[x,y] = Odometer[x,y] + Z_lat[x,y];
            Z_lat[x,y] = 0;

            pop!(w_Bound);

            i1 = min(i1,x-1);
            i2 = max(i2,x+1);
            j1 = min(j1,y-1);
            j2 = max(j2,y+1);

            for i=1:size(v_check,1)
                a = x+v_check[i,1];
                b = y+v_check[i,2];

                if Z_lat[a,b]!=0
                    is_Interior = ((V_sites[a,b + 1]==true)&&(V_sites[a,b - 1]==true)&&(V_sites[a - 1,b]==true)&&(V_sites[a + 1,b]==true));

                    if (is_Interior==false) && (Z_lat[a,b]<=Norm_mass)
                        continue
                    end

                   if is_Interior == true
                       j = FindCoordinate(w_Bound,a,b);

                       if j > 0
                          splice!(w_Bound,j); #removes the element with index j
                      end
                  end

                 # case 2. if (a,b) needs to be walked out,
                 # but it's not in the boundary stack. we check, to see if it needs to be added
                 if is_Interior==false
                    j = FindCoordinate(w_Bound,a,b);
                    if j == -1
                       push!(w_Bound,L_coord(a,b));
                    end
                 end

               end # Z_lat is non zero
             end # end of going over v_check

             top_ = length(w_Bound);

        end # end of while

        t2 = time_ns();
        print( "time elapsed: " , (t2 - t1)/1.0e9, "\n" );

        # =  = = we now move the interior

        println("Moving interiors.....");
        println("..... Step1. Populating the interior....");
        t1 = time_ns();

        # we first populate the interior
        w_Int=L_coord[];

        for i=i1:i2
            for j=j1:j2
                is_Interior = ((V_sites[i,j + 1]==true)&&(V_sites[i,j - 1]==true)&&(V_sites[i - 1,j]==true)&&(V_sites[i + 1,j]==true));
                if is_Interior == true
                   k = FindCoordinate(w_Int,i,j);
                      if k==-1
                          push!(w_Int,L_coord(i,j));
                      end
                end
            end
        end

        t2 = time_ns();
        print( "time elapsed: " , (t2 - t1)/1.0e9, "\n" );
        println("..... Step2. Creating transformations....");

        t1 = time_ns();

        Transforms_ = createTransforms(w_Int); #Transforms_ is a pair of 2 matrices [T_matrix, Delta_matrix] in this order
        # T_matrix is the coeff  matrix of the mass transform for the interior walk
        # Delta_matrix is the coeff matrix of the mass distribution to boundary for the interior walk

        weights_ = zeros(Float64, length(w_Int),1 );
        for i = 1:length(w_Int)
          weights_[i]=Z_lat[w_Int[i].x, w_Int[i].y];
        end

        t2 = time_ns();
        print( "time elapsed: " , (t2 - t1)/1.0e9, "\n" );

        println("..... Step3. Computing distributed weigths with inversion operator....");
        t1 = time_ns();
        id_M = zeros(Float64, length(w_Int), length(w_Int) );
        for ii = 1:length(w_Int)
           id_M[ii,ii] = 1;
        end

        weights_=(Transforms_[2]*(inv( id_M  - Transforms_[1]   )) )*weights_;

        t2 = time_ns();
        print( "time elapsed: " , (t2 - t1)/1.0e9, "\n" );

        println("..... Step4. Distributing the weigths....");

        t1 = time_ns();

        for i = 1:length(w_Int)
            a = w_Int[i].x;
            b = w_Int[i].y;

            Z_lat[a+1,b] = Z_lat[a+1,b] + weights_[i];
            Z_lat[a-1,b] = Z_lat[a-1,b] + weights_[i];
            Z_lat[a,b+1] = Z_lat[a,b+1] + weights_[i];
            Z_lat[a,b-1] = Z_lat[a,b-1] + weights_[i];

            Odometer[a,b] = Odometer[a,b]+4*weights_[i];
         end

         for i=1:length(w_Int)
             a = w_Int[i].x;
             b = w_Int[i].y;
             Z_lat[a,b]=0;
         end

        t2 = time_ns();
        print( "time elapsed: " , (t2 - t1)/1.0e9,"\n" );

        # it is left to re-populate the boundary

        println("..... Step5. Populating the new boundary....");

        t1 = time_ns();

        for i=1:length(w_Int)
           #checking neighbors of interior points
            for k = 1:4
               a = w_Int[i].x + dx[k];
               b = w_Int[i].y + dy[k];

              is_Interior = ((V_sites[a,b + 1]==true)&&(V_sites[a,b - 1]==true)&&(V_sites[a - 1,b] == true)&&(V_sites[a + 1,b]==true));

              if (is_Interior == false)&&(Z_lat[a,b] > Norm_mass)
                  j=FindCoordinate(w_Bound,a,b);
                  if j == -1
                      push!(w_Bound,L_coord(a,b));
                  end
              end
            end
        end # for loop over interior

        t2 = time_ns();
        print( "time elapsed: " , (t2 - t1)/1.0e9, "\n" );

        w_Int = L_coord[];

        if length(w_Bound) == 0
           break;  # WE're done!!!
        end

    end #end of the main while


    Z_lat = output.printIntoFile(Z_lat,  0, string("BSand_Z_", T_mass) );
    Odometer = output.printIntoFile(Odometer, 1,string("BSand_OD_", T_mass)  );

    output.saveAsGrayImage(Z_lat, string("BSand_Z_", T_mass), 20, 2);
    output.saveAsGrayImage(Odometer, string("BSand_OD_", T_mass), 20, 2);

    # for the total elapsed time, it's better to use the @time macros on the main call


    return Z_lat, Odometer; # these are trimmed in output module

end # end of function moveSand




function createGraph(coord_List::Array{L_coord,1})
  # given a list of lattice sites, creates a graph based on the lattice neighborhood
  # e.g. (0,0) will be connected to all of (0,1), (1,0), (-1,0) , (0, -1)
  # the graph representation is of the form i -> array_of_integers,  where i is the index of the site
  # and the array_of_integers is the list of indices of the neighbors of i

  N = length(coord_List);

  res_= Graph_rep[]; # the result

  for i=1:N
    push!(res_, Graph_rep(Int[]));
  end


  for i=1:N
    for j=1:N

      if (i!=j) && (abs(coord_List[i].x-coord_List[j].x) + abs(coord_List[i].y-coord_List[j].y)==1 )
        push!(res_[i].neigh,j);
      end

    end
  end

  return res_;
end

function createTransforms(coord_List::Array{L_coord,1})
  # T_matrix=zeros(Float64,1,1);  # coeff of the mass transform for the interior walk
  # Delta_matrix=zeros(Float64,1,1); # coeff of the mass distribution to boundary for the interior walk
  # definitions of these matrices are explained in Section 3 of our paper, see Algorithm 2

  N=length(coord_List);
  Delta_matrix=zeros(Float64,N,N);

  if N==1
    Delta_matrix[1,1]=0.25;
  end

  G = createGraph(coord_List);

  T_matrix = zeros(Float64, N, N);
  for i = 1:N
    T_matrix[i,i] = 1; # make T_matrix as an identity matrix
  end


  #println("........... Step2b. Creating the matrices....");

  for i=1:N

    for p=1:N
      Delta_matrix[i,p]=Delta_matrix[i,p]+T_matrix[i,p]/4;
    end


    for p=1:N
      for j=1:length(G[i].neigh)
        k=G[i].neigh[j];
        T_matrix[k,p]=T_matrix[k,p]+T_matrix[i,p]/4;
      end
      T_matrix[i,p]=0;

    end

  end

  return T_matrix, Delta_matrix;
end


function moveSand_Standard(T_mass, max_iter_Count = 50000)
    # performes divisible sandpile process in a standard way
    # each time checking for lattice sites which need to move and adding those in the moving stack

    Z_lat = zeros(2*N_size+1,2*N_size+1);   #Array(Float64, 2*N_size+1,2*N_size+1); models the Z^2 lattice
    V_sites = falses(2*N_size+1,2*N_size+1);  # all sites which are visited by the sandpile process, are being marked here
    Odometer = zeros(2*N_size+1, 2*N_size+1);


    w_Int = L_coord[];  # the stack of lattice sites which need to move (topple)

    Norm_mass = T_mass^0.5;  # the maximal allowed mass on the boundary
    Z_lat[N_size+1,N_size+1] = T_mass; # the initial mass is set at the origin
    V_sites[N_size + 1 , N_size + 1 ] = true;
    push!(w_Int,L_coord(N_size+1,N_size+1));

    i, j = 1, 0;
    #typed_Zeros = typeof(Z_lat[1,1])(0);

    i1, i2, j1, j2 = N_size+1, N_size+1, N_size+1, N_size+1;

    iter_ = 0;
    ex_mass_max = 0;
    ex_mass_max_new = 0;
    ex_mass = 0;

    while iter_ < max_iter_Count

        iter_ = iter_ + 1;
        if (iter_%1000 == 0)
            println("Iteration N:: ", iter_);
            println("The boundaries are:: ", (i2-i1+1),"x",(j2-j1+1), "\n");
        end


        ex_mass_max = 0;

        for  i = 1:length(w_Int)
           x, y = w_Int[i].x, w_Int[i].y ;
           ex_mass_max = max(Z_lat[x,y] , ex_mass_max);
        end

      for tt=1:1500
          for i=1:length(w_Int)
               # the toppling
               x, y = w_Int[i].x, w_Int[i].y;
               if Z_lat[x,y] == 0.0
                  continue
               end

               dist_mass = Z_lat[x,y]/4;

               Z_lat[x,y + 1] += dist_mass;
               Z_lat[x,y - 1] += dist_mass;
               Z_lat[x - 1,y] += dist_mass;
               Z_lat[x + 1,y] += dist_mass;

               Odometer[x,y] += Z_lat[x,y];

               Z_lat[x,y] = 0 ;
          end

          ex_mass_max_new = 0 ;
          for i = 1:length(w_Int)
              x, y = w_Int[i].x, w_Int[i].y;
              ex_mass_max_new = max(Z_lat[x,y]  - 1, ex_mass_max_new);
          end
          if ex_mass_max_new < 0.5*ex_mass_max
             break;
          end
       end # interior cycles


        for i = 1:length(w_Int)
            x , y = w_Int[i].x, w_Int[i].y;

            i1, i2, j1, j2 = min(i1,x-1), max(i2,x+1), min(j1,y-1), max(j2,y+1);

            V_sites[x, y + 1] = true;
            V_sites[x, y - 1] = true;
            V_sites[x - 1, y] = true;
            V_sites[x + 1, y] = true;
        end

        w_Int = L_coord[];
        q = false;

        for i=i1:i2
            for j=j1:j2
                is_Interior = ((V_sites[i,j + 1]==true)&&(V_sites[i,j- 1]==true)&&(V_sites[i - 1,j] == true)&&(V_sites[i + 1,j]==true));
                if (is_Interior)
                    push!(w_Int, L_coord(i,j));
                    if Z_lat[i,j] > 0
                      q = true;
                    end
                else
                    if Z_lat[i,j] > Norm_mass
                        push!(w_Int, L_coord(i,j));
                        q = true;
                    end
                end
            end
        end

        if q == false
          break;  # We are done
        end

    end #end of the main while


  println("saving into CSV...")
  Z_lat = output.printIntoFile(Z_lat,  0, string("BSand_Z_", T_mass), true );
  Odometer = output.printIntoFile(Odometer, 1,string("BSand_OD_", T_mass) , true );

  println("saving as Gray images...")
  output.saveAsGrayImage(Z_lat, string("BSand_Z_", T_mass), 20, 0, true);
  output.saveAsGrayImage(Odometer, string("BSand_OD_", T_mass), 20, 0, true);

  return Z_lat, Odometer; # these are trimmed in output module

end # end of function moveSand_Standard


end # end of the module
