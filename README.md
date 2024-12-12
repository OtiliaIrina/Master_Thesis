# Master_Thesis

The notebook class.ipynb is my most developped one, with classes to cresate my own packages.
My trial notebooks are m-vs-z.ipynb (spilts a catalogue into subsamples and comuptes r0), Final-2PCF.ipynb (computes w(theta) by hand and with treecoor, and determines the galaxy bias).
Cosmic variance is in developement to anaylse all the 4 deep fields, to compute the sample variance.


## Theory: 
### Angular Correlation Function (ACF)

The Angular Correlation Function (ACF), a statistical tool used in cosmology to analyze the clustering of galaxies across the universe. It is a way to quantify how likely it is to find two galaxies at a particular angular separation compared to a random distribution. It reveals the non-uniform distribution of galaxies in the universe, helping us understand the large-scale structure and processes at play.

**How Does it Work?**

1. **Galaxy Catalog:** Astronomers create a catalog listing the positions (on the celestial sphere) of observed galaxies.
2. **Pair Counting:** The number of galaxy pairs within specific angular separations is calculated.
3. **Comparison to Random:** This count is compared to the expected number of pairs in a random distribution of galaxies.
4. **Correlation Function (w(θ)):** The ACF quantifies the excess or deficit of galaxy pairs at different angular separations (θ).

**What Can We Learn from the ACF?**

* **Large-Scale Structure:** A high w(θ) at large angular separations indicates galaxies clustered together across vast distances.
* **Small-Scale Structure:** A high w(θ) at small separations suggests galaxies forming groups and clusters.
* **Cosmological Parameters:** By studying the ACF at various scales, scientists can estimate parameters like matter density and dark energy density.


**Applications of the ACF:**

* **Understanding Cosmic Structure Formation:** The ACF provides insights into how galaxies and large-scale structures formed.
* **Testing Cosmological Models:** By comparing the observed ACF with predictions from various models, scientists can validate these models.
* **Dark Matter and Dark Energy:** The ACF helps constrain the properties of dark matter and dark energy, which dominate the universe's mass and energy content.

**r0: A Measure of Clustering Strength**

r0 is a parameter quantifying how strongly galaxies are clustered. It represents the characteristic scale over which this clustering occurs.

**How r0 Changes:**

* **Redshift Dependence:**
    * **Higher Redshift:** At earlier cosmic times (higher redshifts), galaxies are more clustered due to a denser universe and stronger gravitational forces. This translates to a larger r0.
    * **Lower Redshift:** As the universe expands, clustering weakens, leading to a smaller r0.
* **Stellar Mass Dependence:**
    * **Massive Galaxies:** These tend to be more clustered due to deeper gravitational potentials and residing in denser environments (like clusters). Consequently, r0 is typically higher for massive galaxies.
    * **Less Massive Galaxies:** These are less clustered and found in less dense environments, resulting in a smaller r0.

**Understanding the Implications:**

Studying how r0 varies with redshift and stellar mass reveals valuable insights into:

* **Galaxy Formation and Evolution:** A higher r0 suggests galaxies formed and grew more efficiently in that epoch.
* **Nature of Dark Matter and Dark Energy:** The evolution of r0 with redshift can provide clues about the influence of these mysterious components on cosmic structure.

**Visualizing the Characteristic Scale:**

Imagine a cosmic web where galaxies are connected like nodes in a network. r0 represents the average distance between these nodes. In regions with strong clustering (larger r0), the nodes are closer together. Conversely, weaker clustering (smaller r0) signifies more spread-out nodes. 

This document provides a foundational understanding of the ACF and its role in unraveling the mysteries of the cosmos.