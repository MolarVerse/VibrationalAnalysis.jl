"""
	check_unit(unit::String)

Check if the unit is valid.

# Arguments
- `unit::String`: The unit to check.

# Returns
- `function`: The wavenumber function to use.
"""
function check_unit(unit::String)

	# Check if the unit is valid - lowercase
	unit = lowercase(unit)

	# Check if the unit is valid
	if unit == "kcal"
		return wavenumber_kcal
	elseif unit == "hartree"
		return wavenumber_hartree
	elseif unit == "ev"
		return wavenumber_eV
	else
		throw(ArgumentError("Invalid unit. Options are kcal, hartree and eV."))
	end
end