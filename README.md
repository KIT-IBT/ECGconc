# ECGconc
ECGconc offers ECG-based potassium concentration estimation. 

## The repository contents
*Prepare_Data.m*: Template building and extraction of features from a lead reduced signal. An example with data from [1,2] is taken to offer an example for the workflow.

*Apply_Conc_Model.m*: With this script, the parameterized concentration models (from the folder globalModels) can be applied to the found features. The user should keep in mind that the model still needs concentration measurements from a blood test to apply a patient-specific correction.

*Find_Conc_Model_Pat.m*: Code for the generation of the patient-specific model.

*Find_Conc_Model_Global.m*: Code for the generation of the global model.

## Dependencies
ECGdeli: https://github.com/KIT-IBT/ECGdeli

ECGfeat: https://github.com/KIT-IBT/ECGfeat


## Acknowledgements
Special thanks go to Cristiana Corsi and Stefano Severi from the University of Bologna.

## Final Remarks
This is the result of my PhD project at the Institute of Biomedical Engineering at Karlsruhe Insitute of Technology. Further information on the background of this repository will be available when the thesis is published.

## Literature
[1] R. Bousseljot, D. Kreiseler, and A. Schnabel, “Nutzung der EKG-Signaldatenbank CARDIODAT der PTB über das Internet,” Biomedizinische Technik/Biomedical Engineering, pp. 317–318, Jan. 2009.
[2] A. L. Goldberger, L. A. Amaral, L. Glass, J. M. Hausdorff, P. C. Ivanov, R. G. Mark, J. E. Mietus, G. B. Moody, C. K. Peng, and H. E. Stanley, “PhysioBank, PhysioToolkit, and PhysioNet: components of a new research resource for complex physiologic signals.,” Circulation, vol. 101, no. 23, pp. E215–20, Jun. 2000.
