using DelimitedFiles
using LinearAlgebra

c = 2.99792458e10 # in cm/s https://discourse.julialang.org/t/the-speed-of-light-is-an-integer-why-should-we-care/40108

include("dicts__massesdict.jl")

rst = readdlm("start_h2o.rst")
atomNames = rst[:,1][1:end-1]
nAtoms = size(atomNames)[1]

masses = (x->masses[lowercase(x)]).(atomNames)
masses_repeat = repeat(masses, inner=3)
masses_mat = (masses_repeat * masses_repeat').^(1/2)

hess = readdlm("hess.mat")

H = hess' * hess

v, V = eigen(H)

H_sym = V * Diagonal(sqrt.(abs.(v))) * V'

H_sym_w = H_sym ./ masses_mat

### center of mass
coord = rst[:, 4:6][1:end-1, :]
Rcm = sum(masses .* coord, dims = 1) ./ sum(masses)

### Shift coordinates to center of mass frame
R = coord .- Rcm

### Inertia tensor
x = R[:, 1]
y = R[:, 2]
z = R[:, 3]
Ixx = sum(masses .* (y.^2 + z.^2))
Iyy = sum(masses .* (x.^2 + z.^2))
Izz = sum(masses .* (x.^2 + y.^2))
Ixy = - sum(masses .* (x .* y))
Ixz = - sum(masses .* (x .* z))
Iyz = - sum(masses .* (y .* z))
### Create Inertia Tensor
I = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
eigI, X = eigen(I) # Calculate Eigenvalue EigenVector for Inertia Tensor

### Create Translational Vectros
list = [diagm(0 => [i, i, i]) for i in masses]
translation = sqrt.(mapreduce(permutedims, vcat, list))

### Shift coordinates to center of mass frame
# R = coord .- Rcm
### Rotate coordinates into rotational frame (Eckart Frame)
P = R * X
### Rotational Vectors D4-D6
D4 = ((P[:,2] * X[3, :]' .- P[:,3] * X[2,:]') .* sqrt.(masses))[:]
D5 = ((P[:,3] * X[1, :]' .- P[:,1] * X[3,:]').* sqrt.(masses))[:]
D6 = ((P[:,1] * X[2, :]' .- P[:,2] * X[1,:]') .* sqrt.(masses))[:]

# D1-D3 Vectros NORMILIZED
D1 = translation[:,1] / norm(translation[:,1])
D2 = translation[:,2] / norm(translation[:,2])
D3 = translation[:,3] / norm(translation[:,3])
# D4-D6 Vectros NORMILIZED
D4 /= norm(D4)
D5 /= norm(D5)
D6 /= norm(D6)

# Create D Matrix for water
D = zeros((9,6))
# Add D1-D6 to D
D[:, 1] = D1
D[:, 2] = D2
D[:, 3] = D3
D[:, 4] = D4
D[:, 5] = D5
D[:, 6] = D6

### Gram-Schmidt Diagonalization
D_gs = qr(D).Q

### Transform the Hessian to internal coordinates and diagonalize
f_int= D_gs' * H_sym_w * D_gs
eigVal, eigVec = eigen(f_int)
M = Diagonal(1 ./ sqrt.(masses_mat))
Icart = M * D_gs * eigVec

# Normalize Vectors
N = sqrt.(1 ./ sum(Icart.^2, dims=1))
Inormal = Icart .* N

### Conversion kcal Å^-2 g^-1 to s^-2: 4184 * 10^23 (4184 J/kcal * 10^20 Å^2 / m^2 * 1000 g / kg)
eigenValScaled = eigVal * 4184 * 1.0E23
### frequency
ω = sqrt.(eigenValScaled)
### wavenumber
ν_tilde = 1/(2π * c) * ω

#### Reduced masses of the modes
# red_mass = sum(eigVec.^2 .* masses_repeat, dims=1)

red_mass = N.^2

### Force constant of the modes
### Conversion g mol^-1 s^-2 to mdyn Å^-1: / 6.022 / 1E23 (mol) / 1E3 (kg/g) / 1E2 (mdyn/Å / N/m)
k = ω.^2 .* red_mass' / 6.022 / 1E28

### Charge HARD CODED
q = [ -0.65966 , 0.32983 , 0.32983]

intensities = []
### Infrared Intesity
for i in 1:size(Inormal)[1]
    eigenVec = reshape(Inormal[:, i], 3, :)'
    μ1 = sum((R .+ 0.05eigenVec) .* q, dims=1)
    μ2 = sum((R .- 0.05eigenVec) .* q, dims=1)
    Δμ = μ1 - μ2 # norm
    intensity = sum((Δμ ./ 2norm(0.1eigenVec)).^2) # is not right i think should be ./ norm(0.1eigenVec)
    push!(intensities, intensity)
end

#### Print intensities
writedlm(stdout, hcat(ν_tilde, intensities))

### Print all modes
for i in 1:size(Inormal)[1]

    modeIndex = i
    file = open("mode_jl_$modeIndex.xyz", "w")

    mode = reshape(Inormal[:, modeIndex], 3, :)'

    for (i,α) in enumerate(-0.25:0.01:0.25)
        println(file, nAtoms, "\n")

        for i in 1:nAtoms
            println(file, atomNames[i] , "    ",
            join(R[i, :] .+ (α * mode[i,:]), "   " ))
        end
    end

    close(file)

end
