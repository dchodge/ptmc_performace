library('tidyverse')

clean_data <- function(age_c)
{
    age_lim <- c(age_c, 100)
    raw <- read.csv("./rsv_output_samples.csv", skip = 1, stringsAsFactors = FALSE)
    lol <- slice(raw[421,], 1) #Example of an empty row
    age_g <- c(0:11/12, 1:5, 10, 15, 25, 35, 45, 55, 65, 75)
    age_no <- length(age_g)

    pos_empty <- 0
    for(i in 1:nrow(raw))
    {
      if (all.equal(slice(raw[i,],1), lol)==TRUE)
      {
        pos_empty <- append(pos_empty, i)
      }
    }
    pos_empty <- append(pos_empty, nrow(raw))
    age_col <- c()
    for (i in 1:age_no)
    {
      age_col <- append(age_col, rep(age_g[i], pos_empty[i+1] - pos_empty[i]))
    }
    raw$age_group <- age_col

    pos_drop <- c(pos_empty[2:length(pos_empty)-1], pos_empty[2:length(pos_empty)-1]+1, pos_empty[2:length(pos_empty)-1]+2)
    clean <- raw[-pos_drop,]

    clean <- clean %>% mutate(year=substr(D_WeekNo, 0,4), week=substr(D_WeekNo, 5,6)) %>%
      select(Year=year, Week=week,Age.group=age_group,
                              RSV.Negative=RSV.Negative,
                              RSV.Positive=RSV.Positive) %>%
      mutate(RSV.Negative=case_when(RSV.Negative==""~"0",TRUE~RSV.Negative),
                              RSV.Positive=case_when(RSV.Positive==""~"0",TRUE~RSV.Positive)) %>%
      mutate(Age.group=as.double(Age.group),
             RSV.Negative=as.double(RSV.Negative), RSV.Positive=as.double(RSV.Positive))


    clean_edit <- data.frame(Year=character(),
                             Week=character(),
                             Age.group=numeric(),
                             RSV.Negative=numeric(),
                             RSV.Positive=numeric()
    )
    clean_edit_t <- clean_edit
    clean_edit_age <- clean_edit

    #Clean data so that only 52 weeks (26-25) of each year are included
    for(a in 1:age_no)
    {
      clean_age <- clean %>% filter(Age.group==age_g[a])
      for(i in 1:nrow(clean_age))
      {
        if(clean_age[i,2]=="26"&&clean_age[i,1]<"2017")
        {
          clean_edit <- bind_rows(clean_edit, slice(clean_age, i:(i+51)))
        }
      }
    }
    
    clean_edit$Week <- sapply(c(26:52, rep(1:52,6), 1:25), function(x) toString(x)) %>% factor(levels = sapply(1:52, function(x) toString(x)))
    clean_edit$Year <- sapply(c(rep(2010,27), sapply(2011:2016, function(x) rep(x,52)), rep(2017,25)), function(x) toString(x))
    
    
    #Edit the table to account for the grouping
    for (a in 1:(length(age_lim)-1))
    {
      clean_edit_t <- clean_edit %>% filter((Age.group>=age_lim[a])&(Age.group<age_lim[a+1]))
      clean_edit_t <- summarise(group_by(clean_edit_t, Year, Week),
                  RSV.Negative=sum(RSV.Negative),
                  RSV.Positive=sum(RSV.Positive)
        )
      
      clean_edit_t$Age.group <- rep(a, nrow(clean_edit_t))
      clean_edit_age <- bind_rows(clean_edit_age, clean_edit_t)
    }

    s <- c()

    for (a in 1:(length(age_lim)-1))
    {
      tot_pos_y7 <- clean_edit_age %>% filter(Age.group==a) %>% slice((52*6+1):(52*7)) %>% pull(RSV.Positive) %>% sum
      
      for (y in 1:7)
      {
        tot_pos_yx <- clean_edit_age %>% filter(Age.group==a) %>% slice((52*(y-1)+1):(52*y)) %>% pull(RSV.Positive) %>% sum
        s <- append(s, tot_pos_y7/tot_pos_yx)
      }
    }


    s_edit <-c()
    # Make the scalar programme
    for (l in 1:((length(age_lim)-1)*7))
    {
      s_edit <- append(s_edit, clean_edit_age$RSV.Positive[(52*(l-1)+1):(52*l)]*s[l])
    }
    clean_edit_age$RSV.Positive.S <- s_edit

    clean_spread <- clean_edit_age %>% select(Year, Week, Age.group, RSV.Positive.S)# %>% spread(Age.group, RSV.Positive.S, drop = TRUE) #%>% select(sapply(1:25, function(x) toString(x)))
    clean_spread$Week <- sapply(clean_spread$Week, as.numeric) 
    clean_spread <- clean_spread %>% spread(Age.group, RSV.Positive.S, drop = TRUE) %>% select(sapply(1:(length(age_lim) - 1), function(x) toString(x)))
}
