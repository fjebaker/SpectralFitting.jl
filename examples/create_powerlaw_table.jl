"""
Example: Creating a simple power law table model

This script demonstrates how to create an XSPEC-compatible table model
for a power law with a single parameter (photon index).

The power law is defined as: F(E) = E^(-Γ)
where Γ is the photon index parameter.
"""

using SpectralFitting

# Define energy grid: 0.1 to 20 keV with 500 bins (201 edges)
energy_bins = collect(range(0.1, 20.0, length=501))
println("Energy grid: $(length(energy_bins)) edges ($(length(energy_bins)-1) bins)")
println("  Range: $(first(energy_bins)) to $(last(energy_bins)) keV")

# Define photon index grid: 1.5 to 2.5 with 0.1 step size
photon_index_grid = collect(1.5:0.1:2.5)
param_grids = (photon_index_grid,)
println("\nPhoton index grid: $(photon_index_grid)")
println("  Number of parameter values: $(length(photon_index_grid))")

# Calculate spectra for each photon index
n_energy_bins = length(energy_bins) - 1
n_spectra = length(photon_index_grid)
spectra = zeros(Float64, n_energy_bins, n_spectra)

# Compute mid-point energies for evaluation
E_mid = (energy_bins[1:end-1] .+ energy_bins[2:end]) ./ 2

# Fill spectra array: F(E) = E^(-Γ)
for (idx, gamma) in enumerate(photon_index_grid)
    spectra[:, idx] = E_mid .^ (-gamma)
end

println("\nSpectra array size: $(size(spectra))")
println("  ($(size(spectra, 1)) energy bins × $(size(spectra, 2)) parameter combinations)")

# Write to FITS file
output_file = "powerlaw_table.fits"
write_table_model(
    output_file,
    energy_bins,
    param_grids,
    spectra;
    model_name = "POWLAW1D",
    param_names = ["PhoIndex"],
    param_units = [""],
    model_units = "photons/cm^2/s",
    additive = true,
    lowE_units = "keV",
    highE_units = "keV"
)

println("\n✓ Table model written to: $(output_file)")

# Verify the file structure
println("\nVerifying the written file...")
using FITSIO
f_check = FITS(output_file)
println("  Number of HDUs: $(length(f_check))")
println("  HDU 2 (PARAMETERS): $(read_key(f_check[2], "EXTNAME")[1])")
println("  HDU 3 (ENERGIES): $(read_key(f_check[3], "EXTNAME")[1])")
println("  HDU 4 (SPECTRA): $(read_key(f_check[4], "EXTNAME")[1])")

# Check parameter info
param_names_read = read(f_check[2], "NAME")
param_numbvals = read(f_check[2], "NUMBVALS")
println("\n  Parameters:")
for i in 1:length(param_names_read)
    println("    $(param_names_read[i]): $(param_numbvals[i]) values")
end

# Check energy bins
energy_lo_read = read(f_check[3], "ENERG_LO")
energy_hi_read = read(f_check[3], "ENERG_HI")
println("\n  Energy bins: $(length(energy_lo_read))")
println("    Range: $(first(energy_lo_read)) to $(last(energy_hi_read)) keV")

# Check spectra
spectra_read = read(f_check[4], "INTPSPEC")
println("\n  Spectra array: $(size(spectra_read))")

close(f_check)

println("\n✓ File structure verification complete!")
println("\nThe table model file has been successfully created.")
println("You can use this in XSPEC by loading it as a table model.")
