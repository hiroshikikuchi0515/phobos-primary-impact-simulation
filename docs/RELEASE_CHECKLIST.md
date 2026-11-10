# Before making a public Icarus-linked software release

- [ ] Confirm the author(s), institutional requirements, and whether this source may be redistributed publicly.
- [ ] Choose and add an explicit source-code license after rights review; decide separately how derived data are licensed.
- [ ] Audit the entire repository for hard-coded home paths, account details, credentials, private URLs, or other information that must not be published.
- [ ] Confirm exact C++ source hash matches the submitted production run; archive source and compilation options.
- [ ] Confirm which validation and analysis scripts correspond to the final accepted manuscript.
- [ ] Include a run manifest for all final 36 cases and a smaller, practical example that users can execute.
- [ ] Verify the provenance, scientific permissions, and coordinates/units of the shape model and observational catalogs.
- [ ] Provide a clean explanation for obtaining input materials that cannot be redistributed.
- [ ] Provide aggregated outputs and checksums sufficient to regenerate manuscript tables and plots (or clearly state which cannot be regenerated).
- [ ] Verify the code builds on a clean machine and run the smallest documented test.
- [ ] Finalize CITATION.cff, including any collaborators and the version number.
- [ ] Only after checking files: publish the public repository, create a tagged GitHub Release, and if desired archive it on Zenodo.
- [ ] Use the immutable release URL or Zenodo DOI in the paper; do not cite only a moving `main` branch.
- [ ] Update the manuscript Data/Code Availability statement; remove all AUTHOR ACTION placeholders.
