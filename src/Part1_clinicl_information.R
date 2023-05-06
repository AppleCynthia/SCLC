


############## visualization of clinical data
clinic_data <- data.frame(
							Patient_ID = c(paste0('P', seq(1,10))),
							Smoking_year = c(40, 40, 40, 0 ,20,30,0,30,40,0),
							Gender = c('Male', 'Male', 'Male', 'Male',
								'Male', 'Male', 'Male', 'Male', 'Male','Male'),
							Age = c(62, 61, 64,66,56,50,39,65,65,67),
							Treatment = c('Chemotherapy','Treatment_naive', 'Chemotherapy', 
								'Treatment_naive', 'Chemotherapy', 'Treatment_naive', 
								'Chemotherapy', 'Chemotherapy','Treatment_naive','Treatment_naive' ),
							Result = c('Resistance','Treatment_naive', 'No_Responce', 
							           'Treatment_naive', 'Responce', 'Treatment_naive', 
							           'Responce', 'No_Responce','Treatment_naive','Treatment_naive' ),
							Metastasis = c(rep('Metastasis', 8), 'Non_metastasis','Metastasis'),
							Disease_stage = c('IIIA', 'IIIA', 'IIIA', 'IIB', 'IIIA', 'IIIA',
								'IIIB', 'IIIB', 'IIA','IIIA'))

df <- melt(clinic_data, id = 'Patient_ID')

library(ggplot2)
library(ggsci)

result_colors <- c('#FF8000','#FF69B4','#EECD55','#00A087FF')
names(result_colors) <- c('Responce','Resistance','No_Responce','Treatment_naive')

treatment_colors <- c('#3C5488FF','#00A087FF')
names(treatment_colors) <- c('Chemotherapy', 'Treatment_naive')

stage_colors <- rev(pal_locuszoom()(4))
names(stage_colors) <- c('IIA','IIIA','IIB','IIIB')

gender_colors <- c("#40a9bc")
names(gender_colors ) <- 'male'

meta_colors <- c("red","grey") 
names(meta_colors) <- c('Metastasis','Non_metastasis')

ggplot(df, aes(x = Patient_ID, y = variable, label = value)) + geom_tile(aes(fill = value)) + 
geom_text(aes(label=ifelse(value < 70 , as.character(value), '')),
	color ='grey',size=5)+
scale_fill_manual(values= c(rep('white',13), '#3C5488FF',rev(pal_locuszoom()(4)),'#40a9bc','red','#EECD55','gray','#FF69B4','#FF8000', '#00A087FF')) + 

  #scale_fill_manual(values= c(result_colors,treatment_colors,stage_colors,gender_colors,meta_colors))+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0.5, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position='bottom',
                        legend.text=element_text(size=12),
                        axis.title = element_blank(),
                        axis.text.x = element_text(),
                        axis.text = element_text(size=12))+
scale_x_discrete(limits = c(paste0('P', seq(1,10))))

ggsave(file = '/Part1/clinicdata.pdf',
	width = 8, height = 4)

