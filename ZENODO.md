# Publishing (GitHub → Zenodo DOI)

This directory is the **archive tree** used for release: it contains code and
documentation only — no data and no figures.

## 1. Update the GitHub repository

```bash
# Inside the archive directory: initialise / attach the remote (first time only)
cd zenodo_archive/dlbcl-hmgb1-havcr2-mafb-analysis-v1.0.0
git init -b main
git remote add origin https://github.com/tjwu87/dlbcl-hmgb1-havcr2-mafb-analysis.git

git add -A
git commit -m "v1.0.0: reorganise by paper section (Results 1-6); refresh analysis code"
git push -u origin main --force      # overwrite the remote default branch, wiping the old content
```

> If you want the old code to be **completely invisible** on GitHub (including its
> commit history), the cleanest route is to delete the repository, recreate it under
> the same name, and then push this directory.

## 2. Update Zenodo while keeping the same DOI (no new version)

Zenodo only allows self-service file replacement **within 30 days of publication**:

1. Open the record page → **Edit** → expand **Edit files** → **Edit published files**
2. Replace the file (upload the new ZIP, delete the old ZIP)
3. Click **Publish**

Official statement: *"Publishing the draft will not change the DOI."*
(The draft must be published within 45 days of the original publication.)

⚠️ **Do not create another GitHub Release.** This record was created by the GitHub
Release integration; every new Release makes Zenodo mint a **new version + new DOI**.
Updating the repository's default branch does not affect Zenodo.

If "Edit published files" is unavailable (GitHub-linked records can be restricted),
the fallback is: on the record page, set the old file to **restricted access**, so the
DOI record and its public metadata remain but the old code cannot be downloaded.
Contact Zenodo support if necessary.

## 3. If you later need to publish a proper new version

Zenodo record page → **New version** → upload this ZIP → Publish.
This mints a new **Version DOI**; the original **Concept DOI** is unchanged and always
points at the latest version.

## 4. Archive integrity

See `MANIFEST.md` (per-file SHA256).
