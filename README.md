# GeneAnalysis
Documentation 
Project description
 
This graphical user interface (GUI) allows the user to visualise and interact with dynamic graphs that change based on user input. The user can choose different genes and acceptor sites to explore. This Python program uses Plotly Dash to create the interactive graphs with callbacks and triggers, and SQLite to query the database.
 
1.1. Purpose
 
This project investigates the impact of splicing modulator compounds (SMCs) on the proportion of exons expressed in a gene. SMCs target cryptic exons, which are not usually transcribed, and can include them in final RNA transcripts. We are interested in this as the inclusion of a cryptic exon has potential to introduce a premature stop codon into a transcript, triggering nonsense-mediated decay and reducing gene expression. This is particularly important for exploring cancer therapeutics and drug targets. Examples of splicing modulator compounds include Branaplam and Risdiplam.
 
1.2. Goals
 
The goal of this project is to create a user friendly tool that displays plots and visualisations, making the data easy to interpret. 
 
1.3. Intended audience
 
This tool is designed to be used by science lab research and bioinformatics staff that wish to investigate cryptic exons and splicing modulator compounds. 
 
Installation 
2.1. Cloning repository 
 
Cloning the Github repository onto your local computer
(Can refer to: https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
 
  1) On the main page, above the files, click <> Code
  2) Copy the URL
    To clone the repository using HTTPS, under "HTTPS", click the square copy button next to the URL. (https://github.com/Jaime-T/GeneAnalysis.git)
    To clone the repository using an SSH key, including a certificate issued by your organization's SSH certificate authority, click SSH, then click the copy button
    To clone a repository using GitHub CLI, click GitHub CLI, then click the copy button.
  3) Open your terminal
  4) Change directory to where you want the cloned directory to be
  5) Type 'git clone' and then paste the URL you copied in Step 2. Then press Enter to create your local clone
 
2.2. Dependencies
 
The following dependencies are used:  
  pandas
  requests
  numpy
  matplotlib.pyplot
  seaborn
  plotly.express
  dash
  dash_bio
  scipy.special
  sklearn.decomposition.PCA
  plotly.graph_objects
 
These need to be installed separately with a package manager such as pip, e.g. 'pip install pandas'
 
2.4. How to run the program
 
  1) Clone the Github repository (See Section 2.1.)
     
  2) Change the squlite database query path to access the location of where you have the file. This will need to be done in the ‘get_gene_acceptor_data’ and ‘get_raw_expression_data’   function. For example: conn = sqlite3.connect('/Users/folder/example/jc_custom_STARjunc.sqlite')
     
  3)  Run the program by typing into the terminal 'python3 squlite_visualization.py' and press Enter
     
  4) Wait for the program to load. Copy the http URL and paste it into a web browser
 
Key features

3.1. Drop down selection boxes
 
There are multiple drop down boxes in the webapp which allow the user to choose a value.
The first dropdown menu 'Select Gene' allows the user to choose which gene is studied in the needle plot below. The default value is SMN2.
The second dropdown menu 'Show or hide range slider' enables or disables the range slider for the needle plot. The slider is below the needle plot. The range slider allows you to see a window of acceptors.
The third dropdown menu 'maximum needle value' allows the user to pick a value that caps the value of the needles.
The fourth dropdown menu 'Needle value options' allows the user to choose the actual needle value or the logarithm of the needle value, which is useful when the needles are really large, and exceed the maximum needle values.

3.2. Acceptor site Needle Plot
 
The first needleplot displays the exon acceptor sites along the x-axis and the needle value along the y-axis. The needle value represents the significance of the reads at that acceptor. If the treatments have had an impact on the expression, compared to the controls, then the needle value will be higher, allowing the user to identify exons that are targeted by the splicing modulator compounds.
 
The interactive nature of this plot allows the user to click on a needle point (blue dot), which will show a heatmap below with more information about the corresponding acceptor site chosen.
 
For example, the acceptor 70077018 has the highest needle value for gene SMN2. Click on the blue dot to show its heatmap.
 

3.3. Heatmap
 
The heatmap is specific to each gene and acceptor chosen. The heatmap displays the donor sites along the x-axis and the sample description along the y-axis. For each sample row, the proportion of reads across the different donor sites is represented by the colour scale on the right hand side.  The interesting samples have some proportion of reads in the donor site columns that are 'non-canonical', meaning an alternative splicing event is occurring.
 
The information for each cell in the heatmap is displayed when hovering your mouse over the heatmap. The sample label,  donor site, proportion of reads for its row, and raw read count is shown. 
 
 
3.4. Project-based Needle Plot
 
The needle plot is specific to each gene, which is selected at the top of the page. The needle plot has the project ID's along the x-axis and the needle value scores on the y-axis. This feature displays which projects have samples that belong to different Dirichlet distributions, meaning that there is variance in reads between samples in the project, so different splicing events are occurring. The projects with the highest needle values are the ones that are the most affected by the splicing modulator compounds. For example, as shown below, the project SRP334251 has the largest needle value.
 
 
Usage example

(i) This program can be used to see where cryptic exons exist in a gene and which sites can be targeted with splicing modulator compounds to be used as potential drug targets for treating a disease. For example, if you wanted to focus on Huntington's disease, which is caused by mutation in the HTT gene, select HTT from the gene dropdown menu.  The needle plot will show which acceptor sites are most impacted by the experiments.  By clicking on the blue dot of the highest needle value, the heatmap below will update and show a breakdown of the impact of the experiment on each donor site for that acceptor. From here, the heatmap values can be analysed to determine which cryptic exons were expressed, the proportions of the reads, and which splicing modulator compound treatments were effective in causing the splicing change. 
