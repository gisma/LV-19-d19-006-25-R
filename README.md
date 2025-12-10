# LV-19-d19-006-25 · Source Code

## Course Context & Licensing 

![Course Header "Image source: C.Reudenbach"](https://github.com/gisma-courses/LV-19-d19-006-25/blob/main/images/gw-sp.png) 


This course is part of [gisma spatial science ressources](https://gisma-courses.github.io/gc/) of the [Department of Geography](https://www.uni-marburg.de/fb19).

The course content is developed and hosted on Github. 

The responsibility for the content rests with the instructors. Statements, opinions and/or conclusions are the ones from the instructors and do not necessarily reflect the opinion of the representatives of Marburg University.  

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.



## Repository content and purpose

This repository contains the **technical backbone** of the [LV-19-d19-006-25](https://github.com/gisma-courses/LV-19-d19-006-25) course.

It provides a clean, FAIR-oriented structure for all scripts, functions, and processing workflows used throughout the course.

The goal is *not* to offer a single pipeline, but a **collection of reproducible building blocks** that illustrate how modern environmental data workflows are composed:

* **Earth-observation retrieval** (gdalcubes, CDSE, terra)
* **Terrain and GIS processing** (Rsagacmd, link2GI)
* **Statistical and ML components** (tidyverse, caret)
* **Reusable helper functions** for a consistent project layout

Different levels of abstraction are shown side-by-side so users can compare
high-level convenience wrappers with low-level, fully transparent R workflows.

## Structure

```text
renv/    # renv
src/        # scripts &  functions and backend setup
```

## Intended Use

The repository is designed as a **reference implementation** for students:
clone it, inspect the scripts, and use it as a template for your own analyses.
All content is educational; no single workflow is “the” workflow.

## Usage of `renv`

After cloning the repository:

```R
renv::restore()
```

This will automatically install all required R packages in the correct versions as defined in renv.lock.
Your local system libraries do not matter — renv reconstructs the exact project environment for you.

