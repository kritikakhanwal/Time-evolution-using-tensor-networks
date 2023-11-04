using ITensors
using DelimitedFiles

# I: Initialization of the state.

L = 10                         # no of real space sites.
dt = 0.01                          #time step      

s = siteinds("S=1/2",2*L)                # introduce all site indices 
ancilla_sites  = [s[i] for i=2:2:(2*L)]  # an array of indices of ancilla sites
physical_sites = [s[i] for i=1:2:(2*L)]  # an array of indices of physical sites


psi = ITensor[]             # an empty array to store all entangled states
for i=1:2:(2*L)
    ϕ = MPS(s[i:i+1],["+","Up"])  #create an mps state product state |->|0>
    C_not = op("CX",s[i],s[i+1])  # make CNot gate for this two site
    ϕ = apply(C_not,ϕ)# CNOT|->|0> = (|00> + |11>)/sqrt(2)
    append!(psi,ϕ)
end


ψ0 = MPS(psi); # convert the ITensor vector to MPS, each pair is not linked to other pair yet !  
psi = copy(ψ0);  #psi is our state that we will keep evolving! 
#normalize!(psi);

# "Manufacture" all trotter gates needed : 

# Introduce the hamiltonian
physical_gates = []  # an array of strings that contain info of gates needed and physical site index
ancilla_gates  = []  # an array of strings that contain info of gates needed and ancilla site index

# For physical Sites

for j=1:L-1
    s1 = physical_sites[j]
    s2 = physical_sites[j+1]
    hj = op("Sz",s1) * op("Sz",s2) +
         1/2 * op("S+",s1) * op("S-",s2) +
         1/2 * op("S-",s1) * op("S+",s2)
    Gj = exp(-1.0im*dt/2 * hj)
    push!(physical_gates,Gj)
end
    

#For Ancilla sites

for j=1:L-1
    s1 = ancilla_sites[j]
    s2 = ancilla_sites[j+1]
    hj = op("Sz",s1) * op("Sz",s2) +
         1/2 * op("S+",s1) * op("S-",s2) +
         1/2 * op("S-",s1) * op("S+",s2)
    Gj = exp(1.0im * dt/2 * hj)
    push!(ancilla_gates ,Gj)
    
end

steps = 10
cutoff_svd=1E-8
Mdim = 50
psi = copy(ψ0)


for i=1:steps
    
    #---------------Forward_trotter_Physical_indices----------------------------------

    for i=1:L-1              #physical index 
        j = 2*i              #ancilla index   

        orthogonalize!(psi,j); # shift the orthogonal center at ancilla site

        #The swap of physical and ancilla sites so that two physical sites come next to each other
        # This is done by combining the ancilla and the next physical site first and then doing the svd
        # aissgn the physical index and left link to the left tensor so that physical tensor are next to eachother

        U,S,V = svd(psi[j]*psi[j+1],(siteind(psi,j+1),linkind(psi,j-1));cutoff=cutoff_svd , maxdim=Mdim)
        psi[j] = U*S                        # swap the physical index with ancilla
        psi[j+1] = V

        # Appy trotter gate to to i and i+1 physical site which are next to each other

        W = noprime(physical_gates[i]*contract(psi[j-1]*psi[j]))  


        # after the application of trotter, we change the tensor back to MPS state with two sites using SVD

        if(i==1) 
            U,S,V = svd(W,siteind(psi,j-1);cutoff=cutoff_svd,maxdim=Mdim)
        else
            U,S,V = svd(W,(siteind(psi,j-1),linkind(psi,j-2));cutoff=cutoff_svd,maxdim=Mdim)
        end

        #SVD of the two MPS after trotter gate: 
        # keep the physical index to left tensor and ancilla to right 
        psi[j-1] = U
        psi[j]   = S*V


        # Put the ancilla site back to even index place :
        # This is done by combining the tensor at even position (the right tensor in above process)
        # combine it with ancilla next to it and then swap the indices

        U,S,V = svd(psi[j]*psi[j+1],(siteind(psi,j+1),linkind(psi,j-1));cutoff=cutoff_svd,maxdim=Mdim)
        psi[j] = U 
        psi[j+1] = S*V

    end

    #------------------------Backward_trotter_Physical_indices----------------------------------

    for i=reverse(1:L-1)             #physical index 

        j = 2*i                      #ancilla index   
        orthogonalize!(psi,j); 

        #The swap of physical and ancilla sites
        U,S,V = svd(psi[j]*psi[j+1],(siteind(psi,j+1),linkind(psi,j-1));cutoff=cutoff_svd, maxdim=Mdim)
        psi[j] = U*S                        # swap the physical index with ancilla
        psi[j+1] = V 


        # Appy trotter gate:
        W = noprime(physical_gates[i]*contract(psi[j-1]*psi[j]))         

        if(i==1)
            U,S,V = svd(W,siteind(psi,j-1);cutoff=cutoff_svd,maxdim=Mdim)
        else
            U,S,V = svd(W,(siteind(psi,j-1),linkind(psi,j-2));cutoff=cutoff_svd,maxdim=Mdim)
        end

        #SVD of the two MPS after trotter gate 
        psi[j-1] = U
        psi[j]   = S*V

        #print(psi)
        # Put the ancilla site back to even place :
        U,S,V = svd(psi[j]*psi[j+1],(siteind(psi,j+1),linkind(psi,j-1));cutoff=cutoff_svd, maxdim=Mdim)
        psi[j] = U 
        psi[j+1] = S*V

        #print(psi)
    end

    # We repeat the process as for physical site evolution but just take not of the sites:

    #------------------------Forward_trotter_Ancilla_indices----------------------------------

    for i=1:L-1               #ancilla index 

        j = 2*i+1             # physical index   

        orthogonalize!(psi,j);

        #print(psi)

        #The swap of physical and ancilla sites:
        #here the ancilla indices will be next to each other
        U,S,V = svd(psi[j]*psi[j+1],(siteind(psi,j+1),linkind(psi,j-1));cutoff=cutoff_svd, maxdim=Mdim)
        psi[j] = U*S                        # swap the physical index with ancilla
        psi[j+1] = V 

        #print(psi)
        # Appy trotter gate to the pair of ancilla indices
        W = noprime(ancilla_gates[i]*contract(psi[j-1]*psi[j]))         
        U,S,V = svd(W,(siteind(psi,j-1),linkind(psi,j-2));cutoff=cutoff_svd, maxdim=Mdim)

        #SVD of the two MPS after trotter gate is applied
        psi[j-1] = U
        psi[j]   = S*V

        #print(psi)
        # Put the ancilla site back to even place :
        U,S,V = svd(psi[j]*psi[j+1],(siteind(psi,j+1),linkind(psi,j-1));cutoff=cutoff_svd, maxdim=Mdim);
        psi[j] = U 
        psi[j+1] = S*V

        #print(psi)
    end

    #------------------------Backward_trotter_Ancilla_indices----------------------------------

    for i=reverse(1:L-1)      #ancilla index        

        j = 2*i+1             #physical index
        orthogonalize!(psi,j); 
        #print(psi)
        #The swap of physical and ancilla sites
        U,S,V = svd(psi[j]*psi[j+1],(siteind(psi,j+1),linkind(psi,j-1));cutoff=cutoff_svd, maxdim=Mdim);
        psi[j] = U*S                        # swap the physical index with ancilla
        psi[j+1] = V 

        #print(psi)
        # Appy trotter gate:
        W = noprime(ancilla_gates[i]*contract(psi[j-1]*psi[j]))         

        U,S,V = svd(W,(siteind(psi,j-1),linkind(psi,j-2));cutoff=cutoff_svd, maxdim=Mdim)

        #SVD of the two MPS after trotter gate 
        psi[j-1] = U
        psi[j]   = S*V

        #print(psi)
        # Put the ancilla site back to even place :
        U,S,V = svd(psi[j]*psi[j+1],(siteind(psi,j+1),linkind(psi,j-1));cutoff=cutoff_svd, maxdim=Mdim);
        psi[j] = U 
        psi[j+1] = S*V

        #print(psi)


    end

end







