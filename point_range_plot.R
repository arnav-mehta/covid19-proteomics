# takes a dataframe and a protein and generates point range plots

point_range_plot <- function(df, protein, grouping="WHO.max"){

  if (grouping=="WHO.max"){  
    # group by WHO.max and find summary stats
    data <- subset(df, OlinkID==protein)
    plotdata <- data %>%
      group_by(WHO.max,day) %>%
      summarise(n = n(),
                mean = mean(NPX),
                sd = sd(NPX),
                se = sd / sqrt(n))
    assay <-data$Assay[1]
    
    # plot timecourse by WHO.max
    a<-ggplot(plotdata, aes(x=day, y=mean, group=WHO.max, color=WHO.max)) +
      geom_line(size=1) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.1) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(y="NPX") +
      ggtitle(assay)
  }

  else if (grouping=="severity"){
    # group by severity and find summary stats
    data <- subset(df, OlinkID==protein)
    plotdata <- data %>%
      group_by(severity,day) %>%
      summarise(n = n(),
                mean = mean(NPX),
                sd = sd(NPX),
                se = sd / sqrt(n))
    assay <-data$Assay[1]
    
    # plot timecourse by WHO.max
    a<-ggplot(plotdata, aes(x=day, y=mean, group=severity, color=severity)) +
      geom_line(size=1) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.1) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(y="NPX") +
      ggtitle(assay)
  }
  
  else if (grouping=="age"){
    # group by severity and find summary stats
    data <- subset(df, OlinkID==protein)
    plotdata <- data %>%
      group_by(age.group,day) %>%
      summarise(n = n(),
                mean = mean(NPX),
                sd = sd(NPX),
                se = sd / sqrt(n))
    assay <-data$Assay[1]
    
    # plot timecourse by WHO.max
    a<-ggplot(plotdata, aes(x=day, y=mean, group=age.group, color=age.group)) +
      geom_line(size=1) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.1) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(y="NPX") +
      ggtitle(assay)
  }
  
  else if (grouping=="covid"){
    # group by severity and find summary stats
    data <- subset(df, OlinkID==protein)
    plotdata <- data %>%
      group_by(COVID,day) %>%
      summarise(n = n(),
                mean = mean(NPX),
                sd = sd(NPX),
                se = sd / sqrt(n))
    assay <-data$Assay[1]
    
    # plot timecourse by WHO.max
    a<-ggplot(plotdata, aes(x=day, y=mean, group=COVID, color=COVID)) +
      geom_line(size=1) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.1) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(y="NPX") +
      ggtitle(assay)
  }

  else if (grouping=="Vlcat1"){
    # group by severity and find summary stats
    data <- subset(df, OlinkID==protein)
    plotdata <- data %>%
      group_by(Vlcat1, day, severity) %>%
      summarise(n = n(),
                mean = mean(NPX),
                sd = sd(NPX),
                se = sd / sqrt(n))
    assay <-data$Assay[1]

    # plot timecourse by WHO.max
    a<-ggplot(plotdata, aes(x=day, y=mean, group=Vlcat1, color=Vlcat1)) +
      geom_line(size=1) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.1) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(y="NPX") + facet_wrap(~severity) +
      ggtitle(paste(assay, protein))
  }
  
  else if (grouping=="sev_int"){
    
    # group by severity and find summary stats
    data <- subset(df, OlinkID==protein)
    plotdata <- data %>%
      group_by(sev_int,day) %>%
      summarise(n = n(),
                mean = mean(NPX),
                sd = sd(NPX),
                se = sd / sqrt(n))
    assay <-data$Assay[1]
    
    # plot timecourse by WHO.max breaking patients in WHO 1 into intubated or not
    a<-ggplot(plotdata, aes(x=day, y=mean, group=sev_int, color=sev_int)) +
      geom_line(size=1) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.1) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(y="NPX") +
      ggtitle(paste(assay, protein))
  }
  
  else {print("What grouping would you like?")}
  
  a
  
}
