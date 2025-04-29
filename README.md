# Gully Erosion Susceptibility Modeling with R 

### Authors  
### LINK OF THE PAPER - REFERENCE 
**Christian Conoscenti**, **Grazia Azzara**, **Aleksey Y. Sheshukov**
---

### Description  
This repository contains the script for pixel-scale gully erosion susceptibility modeling. It is designed to be reproducible and adaptable to any region where a reliable gully inventory and an accurate DEM are available.

The script automates the snapping of gully trajectories to DEM-derived flow lines using two user-defined parameters: maximum snap distance (SD) and contributing area threshold (CA). The user can set three different values for each parameter (SD1, SD2, SD3 and CA1, CA2, CA3) and adapt the approach to the specific characteristics of their study area.


The modeling framework was implemented in **R**, an open-source programming language, and uses geospatial and statistical analysis packages.  

### **Dependencies**  
The script was developed using R (R Core Team, 2023). To run the script, you need:  
- **R Software**: Download from [CRAN](https://cran.r-project.org/bin/windows/base/).  
- **R Packages**:  
  - **Geospatial Analysis**:  
    - [WhiteboxTools](https://www.whiteboxgeo.com/) (Wu and Brown, 2022; Lindsay, 2016)  
    - [terra](https://rspatial.org/terra/) (Hijmans, 2023)  
  - **Statistical Analysis**:  
    - [data.table](https://cran.r-project.org/package=data.table) (Barrett et al., 2023; Olsen, 2023)  
    - [splitstackshape](https://cran.r-project.org/package=splitstackshape) (Mahto, 2019)  
    - [car](https://cran.r-project.org/package=car) (Fox and Weisberg, 2019)  
    - [earth](https://cran.r-project.org/package=earth) (Milborrow, 2023)  
    - [caret](https://topepo.github.io/caret/) (Kuhn and Max, 2008)  
    - [pROC](https://cran.r-project.org/package=pROC) (Robin et al., 2011)  

### **Input Data Requirements**  
The script requires the following inputs:  
1. **Study Area**: Shapefile (.shp) of the study area, subdivided into 20 spatially tiles. Each tile is a portion of the study area where gullies have been mapped.
2. **DEM (Digital Elevation Model)**: Raster file of the study area.  
3. **Gully Inventory**: Shapefile containing gully trajectory data.  

### **Usage**  
1. Clone or download this repository.  
2. Ensure all required R packages are installed.  
3. Provide the input files (area boundary, DEM, and gully inventory shapefiles).  
4. Run the script in R to perform the modeling and analysis.  

### **License**  
This project is distributed under the **GNU General Public License v3.0**. 
See the full license [here](https://www.gnu.org/licenses/gpl-3.0.html).  

### **Disclaimer**  
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NONINFRINGEMENT.  
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  

### **Contact**  
Department of Earth and Marine Sciences,  
University of Palermo,  
Via Archirafi 22, Palermo 90123, Italy  

**Email**: christian.conoscenti@unipa.it  

---  
For more details, visit the [GitHub repository](https://github.com/graziaazzara/gully_erosion).

