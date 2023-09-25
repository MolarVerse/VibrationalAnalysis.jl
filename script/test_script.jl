using DelimitedFiles
using LinearAlgebra

c = 2.99792458e10 # in cm/s https://discourse.julialang.org/t/the-speed-of-light-is-an-integer-why-should-we-care/40108

include("dicts__massesdict.jl")

rst = readdlm("ceo2_min.rst")
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

### Inertia tensor
x = coord[:, 1]
y = coord[:, 2]
z = coord[:, 3]
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
R = coord .- Rcm
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
D = zeros((size(D1)[1],6))
# Add D1-D6 to D
D[:, 1] = D1
D[:, 2] = D2
D[:, 3] = D3
D[:, 4] = D4
D[:, 5] = D5
D[:, 6] = D6

### Gram-Schmidt Diagonalization
D_gs = qr(D).Q

### Transform the Hessian to interaln coordinates and diagonalize
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
ω = sqrt.(complex(eigenValScaled))
### wavenumber
ν_tilde = 1/(2π * c) * ω

#### Reduced masses of the modes
# red_mass = sum(eigVec.^2 .* masses_repeat, dims=1)
red_mass = N.^2

# Print reduced masses
# println("Reduced Masses")
# writedlm(stdout, red_mass')

### Force constant of the modes
### Conversion g mol^-1 s^-2 to mdyn Å^-1: / 6.022 / 1E23 (mol) / 1E3 (kg/g) / 1E2 (mdyn/Å / N/m)
k = ω.^2 .* red_mass' / 6.022 / 1E28
# println("Force Constants")
# writedlm(stdout, k)

### Charge HARD CODED ATTENTION!!!!!!! 
charge = Dict("ce" => 1.796238775774, "o" => -0.898119387887)
q = (x->charge[lowercase(x)]).(atomNames)

intensities = []
### Infrared Intesity
for i in 1:size(Inormal)[1]
    eigenVec = reshape(Inormal[:, i], 3, :)'
    # http://thiele.ruc.dk/~spanget/help/g09/k_constants.html
    intensity = sum((sum(eigenVec .* q, dims=1) / 0.2081943  / norm(eigenVec)).^2 / red_mass[i] * 42.2561)
    push!(intensities, intensity)
end

#### Print intensities
println("Intensities")
writedlm(stdout, hcat(ν_tilde[7:end], intensities[7:end]))


# Print all modes
# for i in 1:size(Inormal)[1]

#     modeIndex = i
#     file = open("mode_$modeIndex.xyz", "w")
#     # mode_file = open("vec_$modeIndex.txt", "w")

#     mode = reshape(Inormal[:, modeIndex], 3, :)'
#     # writedlm(mode_file, mode)
#     for (i,α) in enumerate(-1:0.01:1)
#         println(file, nAtoms, "\n")

#         for i in 1:nAtoms
#             println(file, atomNames[i] , "    ",
#             join(R[i, :] .+ (α * mode[i,:]), "   " ))
#         end
#     end
#     # close(mode_file)
#     close(file)

# end

# # Calculate dipole moment from e Å to Debye
# μ = sum(R .* q, dims=1) / 0.2081943
# println("Dipole moment in Debye: ", μ )

# # Calculate quadrupole moment from e Å^2 to Debye Å
# gg = R .* q
# quad = gg' * R / 0.2081943

# println("Quadrupole moment: ", quad)

using Plots
using Trapz
using LaTeXStrings

# add a broadening function
function add_broadening(list_ex_energy, list_osci_strength, line_profile::String="Lorentzian", line_param::Float64=10.0, step::Float64=0.1) where T<:Real
    x_min = 0
    x_max = 1250
    x = x_min:step:x_max
    y = zeros(length(x))

    for xp in 1:length(x)
        for (e, f) in zip(list_ex_energy, list_osci_strength)
            if line_profile == "Gaussian"
                y[xp] += f * exp(-((e - x[xp]) / line_param)^2)
            elseif line_profile == "Lorentzian"
                y[xp] += 0.5 * line_param * f / (pi * ((x[xp] - e)^2 + 0.25 * line_param^2))
            end
        end
    end

    # # normalize the lineshape so that its integral is 1
    # integral = trapz(y, x)
    # y /= integral

    # # scale the lineshape to retain the same peak height
    # max_y = maximum(y)
    # y *= maximum(list_osci_strength) / max_y
    
    return x, y
end

# plot the spectrum
function plotSpectrum(ν, I)
    # scatter(ν, I, label="Spectrum")
    x,y = add_broadening(ν, I)
    plot!(x, y, label="Spectrum", xlabel=L"wavenumber in cm${^-1}$", ylabel=L"Intesity in km mol$^{-1}$")
end

file = open("spectrum_ceo2.txt", "w")

# write the spectrum to a file
writedlm(file, hcat(add_broadening(real.(ν_tilde), intensities)...))

# close the file
close(file)
