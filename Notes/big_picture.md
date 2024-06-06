# Big picture 6.6, beginning of project

Forked and prepared folder structure for developing batchCorr for Bioconductor. Made a tentative vignette which illustrates the big picture and helps us in staying on the same page. Vilhelm thought it useful to understand the context to be able to dig into a project, even if the work itself is fairly mechanical. So Vilhelm prepared a vignette to help us with decisions and grasping the task at hand. Vilhelm proposes a prioritized to-do list, work on which is included in the development branch of the fork. Aim to work narrowly on the tasks, without changing function names and folder names, for example. You can get inspiration from the mia and notame packages, for example. 

There are also some guidelines not mentioned in the tests, tasks for which are included in the to-do list. Vilhelm runs the checks and ensures Bioconductor-compatibility of tasks done so everyone doesn't have to set up the development environment (supposedly works best on Linux). Moreover, many tasks are so small that it's not worth keeping track of who's doing what; Vilhelm can do them.

We'll be forking Vilhelm's repository to the Chalmer's metabolomics GitHub? and using a Trello board, where Vilhelm posts tasks and the files it is concerned with. The tasks are developed in their own branch, and merged with the main branch by Vilhelm to ensure Bioconductor-compatibility.

What:
- present the vignette and your co-operation workflow
- begin by writing unit tests for exported functions 

Why:
- to be able to develop without accidentally breaking anything when refactoring and modifying the code

How:
- agree on the overarching structure of the package, which functions are exported

## Help with development environment/running checks

The following is concerned with software packages. Other types of packages include experiment data, annotation and workflow packages. Much of the content herein is from https://contributions.bioconductor.org/index.html. In brief, the package must pass R CMD build and R CMD check using a recent R-devel. The package must also pass "BiocCheck::BiocCheckGitClone()" and "BiocCheck::BiocCheck('new-package'=TRUE). Address all errors, warnings and notes arising during the check; if they are not resolved, a justification is expected by the Bioconductor team. Prioritize errors over warnings over notes and strategize on the order in which problems are addressed. See the check-specific notes and Git commits for notame (https://github.com/vsuksi/notame) for finer granularity and inspiration on how to fix specific problems. Moreover, there are requirements and guidelines glimpsed from the Bioconductor site above, which are considered manually. A few notes:

- As you modify documentation, you need to run devtools::document() in the project directory to render the corresponding .Rd files
- Run devtools::load_all in the project directory to load the package, including non-exported functions, for development
- Run testthat::test_local(path="tests", reporter=testthat::SummaryReporter$new()) in the project directory to specifically check unit tests and get informative error messages
- Run devtools:run_examples() in the project directory to run examples
- many devtools and testthat have useful arguments to, for example, narrow the scope of the check so it passes quicker
- Log files are available for the R CMD check and BiocCheck::BiocCheck() in the directory from which you launched the checks
- A personal access token needs to be used as as a password when setting up Git for pushing to the repository (at least on GitLab). The personal access token is shown only once when it's created, save it somewhere!

### Instructions for running tests:
1. On the command line, navigate to directory containing the batchCorr directory and run "R CMD build x", where x is the batchCorr directory
- the resulting tar.gz directory also contains the rendered vignette in inst/doc
2. On the command line, navigate to directory containing the built tar.gz file and run "R CMD check x" where x is the built tar.gz file
- Bioconductor uses a specific set of R CMD check flags, which instructions online
3. Launch R in the batchCorr directory, install BiocCheck and run BiocCheck::BiocCheck('new-package' = TRUE)
4. Launch R in the batchCorr directory and run BiocCheck::BiocCheckGitClone()


### Development environment:
Use the development version of Bioconductor. Use the development version of R from mid-October to the Bioconductor release in mid-April, use the latest R release from mid-April to the Bioconductor release in mid-October.

I decided to use Linux (Mint Mint 21.3 Cinnamon) for development; apparently there can be unexpected problems with Windows, for example that the checks which must be passed for inclusion in Bioconductor don't notify of all errors, warnings and notes.

To launch R, do `bash R-4.4.0` or `bash R-devel`.

Set in R_profile for semi-transparency so that, for example, labels are displayed in plots:
```
setHook(packageEvent("grDevices", "onLoad"),
function(...) grDevices::X11.options(type='cairo'))
options(device='x11')
```

Download the Bioconductor-specific R CMD check flags (check.Renviron) from Bioconductor website (https://contributions.bioconductor.org/general.html?q=flags#r-cmd-check-environment) and include it in your home/R directory (instructions on website). Reference this in your .Rprofile file by adding "Sys.setenv(R_CHECK_ENVIRON = "~/R/check.Renviron")" . One could check that it works by modifying the flags and seeing the output. Then add "Sys.setenv(R_CHECK_ENVIRON = "~/R/check.Renviron")" to .Rprofile (probably in the home directory).

I used Git version control, check Git log "GitHub link" for workflow and inspiration for problems. Remember to use the login shell so the library paths are set correctly. There's a GUI option for always using the login shell in Linux Mint.

RTools is not required nor available on Linux and Mac. You may need to install the separate program on Windows.

Other things I've done to have it work:
- Added pdflatex to path
- installed some Linux packages to be able to install R packages

### Install the latest R release version

I followed: https://linuxhint.com/install-r-and-rstudio-linux-mint/. In brief:

```
sudo apt update
sudo apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common -y
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
sudo apt update
sudo apt install r-base -y
sudo -i R
```

- I renamed the R executable (in root/bin) to R 4.3.3 to avoid potential confusion

- Update the R installation: `sudo apt upgrade r-base r-base-dev`

- If you get error "Installation paths not writable" (at the beginning of the output) when installing a package, this is because R was installed with some core packages using admin rights. To circumvent this issue, you could sudo launch R for installing the package, or alternatively, set up R as per the following link: https://stackoverflow.com/questions/41839214/installation-path-not-writable-r-unable-to-update-packages

- Some packages may require that additional Linux packages are installed

### Install the latest R development version
I followed https://www.r-bloggers.com/2015/10/installing-r-devel-on-linux-ubuntu-mint/. In brief:
```
sudo apt-get install subversion ccache
sudo apt-get install xorg-dev
cd ~
mkdir svn
cd svn
svn co https://svn.r-project.org/R/trunk r-devel/R
cd r-devel/R/
mkdir ~/svn/R-devel-build
```

Next, include the build-R-devel script (Supplementary files or link) in the R-devel-build directory (add --with-cairo and --with-libtiff). Then: `sudo bash build-R-devel`. I had to additionally install libcurl for the installation to succeed: `apt install libcurl4-openssl-dev`. Next, include the R-devel script under `/usr/local/bin`.

The R-devel script can also be used to to build and check the package, as we will see later. The installation can be updated to the latest version as follows (with update of packages):
```
cd ~/svn/r-devel/R
svn update
sudo bash build-R-devel
bash R-devel
update.packages(ask=FALSE, checkBuilt = TRUE)
```


