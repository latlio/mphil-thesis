<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://www.github.com/latlio/mphil-thesis">
    <img src= "cambridge-logo.png" width="80" height="100">
  </a>

  <h3 align="center">MPhil Thesis 2020 </h3>

  <p align="center">
    The Association Between a Coronary Artery Disease Polygenic Risk Score and Cardiovascular Disease in Women with Breast Cancer
  </p>
</p>

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
* [Running the Code](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
* [Contributing](#contributing)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)

<!-- ABOUT THE PROJECT -->
## About The Project
This is the accompanying R code for my Master's dissertation at the University of Cambridge. For some brief context, advancements in cancer therapeutics have resulted in increases in cancer-related survival and thus, there is an increased burden of cardiovascular disease in breast cancer survivors. This reflects a growing clinical dilemma of balancing the survival benefits and future cardiotoxic harms of oncotherapies. The study is the first to assess the association between a coronary artery disease-specific polygenic risk score (GRS49K) and incident coronary artery events and risk-stratify female breast cancer survivors. 

For confidentiality purposes, I haven't provided the underlying data; however, the code is made available for those who want to replicate my analyses in the future. I primarily use survival analysis (left-truncated Cox PH regression) as well as evaluate C-indices and net reclassification indices. 

<!-- GETTING STARTED -->
## Getting Started
The entire project is split across two scripts: one to prepare the data and the other to run the analyses. I may change this format in the future and convert the code to <code>drake</code> format.

### Prerequisites
You want to make sure that all your file paths are appropriately renamed and that your R is updated the latest version. 

### Installation
 
1. Clone the repo
```sh
git clone https://github.com/latlio/mphil-thesis.git
```
2. Install necessary R packages
```sh
packages_to_install <- c("tidyverse", "survival", "survminer", "compareGroups", "sjlabelled", "rlang", "naniar", "patchwork",
"corrr", "cmprsk")
install.packages(packages_to_install)
```

<!-- USAGE EXAMPLES -->
## Usage
The analysis is designed to be run using only the DoThesisAnalysis.R script. The first time you run the script, you should specify <code>DoAllAnalysis(F)</code> in order to generate the data and save it to your desired directory. Any subsequent runs, you can simply run <code>DoAllAnalysis()</code> because it will load the data. Because everything is wrapped in functions, I find it personally helpful to add a <code>browser()</code> statement within DoAllAnalysis so that I can work within the function environment and run separate lines of code as needed. 

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<!-- CONTACT -->
## Contact

Email: lathanliu21@gmail.com

Project Link: [https://github.com/latlio/mphil-thesis](https://github.com/latlio/mphil-thesis)

