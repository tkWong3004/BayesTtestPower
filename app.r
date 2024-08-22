library(rootSolve)
library(shiny)
library(Rcpp)
library(BH)

sourceCpp("boost_noncentralt.cpp", cacheDir = "tmp_cache")
sourceCpp("pt.cpp", cacheDir = "tmp_cache")
source("onesample.r")
source("twosample.r")

#########################################################

## shiny app
ui <- 
  navbarPage( id = "id",
              "Bayes Factor Design Analysis for Bayesian t-tests",
              

              tabPanel("One-sample t-test and/or paired t-test",
                       sidebarLayout(sidebarPanel(
                         radioButtons("mode", label = p(strong("BFDA:")), choices = c("for sample size determination","for a fixed sample size"), inline = TRUE),
                         
                         
                         p(strong("Specify Alternative hypothese")),
                         radioButtons("rb", label = h4(HTML(paste0("H", tags$sub("1"), " :"))), choices = c("δ ≠ 0","δ > 0", "δ < 0"), inline = TRUE),
                         
                         
                         
                         conditionalPanel(
                           condition = "input.rb == 'δ = δ₁'",
                           numericInput("h1", label = h4("δ ="), value= 1)),
                         conditionalPanel(
                           condition = "input.rb == 'δ ≠ 0'||input.rb == 'δ > 0'||input.rb == 'δ < 0'",
                           radioButtons("model", label = p(strong("Specify analysis prior for δ under H₁:")), choices = c("Cauchy","Normal", "t-student","Non-local"), inline = TRUE)),
                           
                           conditionalPanel(
                             condition = "(input.model == 'Cauchy') && (input.rb == 'δ ≠ 0'||input.rb == 'δ > 0'||input.rb == 'δ < 0')",
                             div(style="display: inline-block; width: 100px;",
                                 sliderInput("l_c", "location:", value = 0,min=-2,max=2,ticks = FALSE,step = .01)),
                             div(style="display: inline-block; width: 100px;",
                                 sliderInput("s_c", "scale:", value = .707,min=.1,max=3,ticks = FALSE,step = .001)),
                             
                           ),
                         conditionalPanel(
                           condition = "(input.model == 'Normal') && (input.rb == 'δ ≠ 0'||input.rb == 'δ > 0'||input.rb == 'δ < 0')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_n", "Mean:", value = 0,max=2,min=-2,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("s_n", "SD:", value = 1,max=3,min=.1,ticks = FALSE,step = .001)),
                           
                         ),
                         conditionalPanel(
                           condition = "(input.model == 't-student') && (input.rb == 'δ ≠ 0'||input.rb == 'δ > 0'||input.rb == 'δ < 0')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_t", "location:", value = 0,max= 2,min=-2,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("s_t", "scale: ", value = .707,min=.1,max=2,ticks = FALSE,step = .001)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("df_t", "df: ", value = 1,min=1,max = 50,ticks = FALSE,step = 1))
                           
                         ),
                         conditionalPanel(
                           condition = "(input.model == 'Non-local') && (input.rb == 'δ ≠ 0'||input.rb == 'δ > 0'||input.rb == 'δ < 0')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_nlp", "location:", value = 0,max=2,min=-2,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 150px;",
                               sliderInput("s_nlp", "scale (mode = location ± scale * √2 ): ", value = .045,max=1,min=.045,ticks = FALSE,step = .001))
                         ),
                         
                         sliderInput("de", p(strong("Specify the bound of compelling evidence BF:")), min=1,max = 20,value = 3),
                        
                         conditionalPanel(
                           condition = "input.mode == 'for sample size determination'",
                           sliderInput("pro", p(strong("Specify the desired probability of true positive evidence")),
                                       min =.5, max = .99, value = .8,step = .01)),
                         
                         conditionalPanel(
                           condition = "input.mode == 'for a fixed sample size'",
                           numericInput("N", p(strong("Sample Size")),
                                        value = 2, min = 2, max = 100000, step = 1))
                         
                     
                         ,
                         
                         
                         radioButtons("daa", label = p(strong("Design prior is the same as analysis prior :")), choices = c("Yes","No"), inline = TRUE),
                         
                         conditionalPanel(
                           condition = "(input.daa == 'No')",
                           radioButtons("model_daa", label = p(strong("Specify the model for design prior:")), choices = c("Cauchy","Normal", "t-student","Non-local"), inline = TRUE)),
                         
                         conditionalPanel(
                           condition = "(input.model_daa == 'Cauchy') &&  (input.daa == 'No')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_c_daa", "location:", value = 0,min=-2,max=2,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("s_c_daa", "scale:", value = .707,min=.1,max=3,ticks = FALSE,step = .001)),
                           
                         ),
                         conditionalPanel(
                           condition = "(input.model_daa == 'Normal') && (input.daa == 'No')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_n_daa", "Mean:", value = 0,max=2,min=-2,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("s_n_daa", "SD:", value = 1,max=3,min=.1,ticks = FALSE,step = .001)),
                           
                         ),
                         conditionalPanel(
                           condition = "(input.model_daa == 't-student') && (input.daa == 'No')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_t_daa", "location:", value = 0,max= 2,min=-2,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("s_t_daa", "scale: ", value = .707,min=.1,max=2,ticks = FALSE,step = .001)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("df_t_daa", "df: ", value = 1,min=1,max = 50,ticks = FALSE,step = 1))
                           
                         ),
                         conditionalPanel(
                           condition = "(input.model_daa == 'Non-local') && (input.daa == 'No')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_nlp_daa", "location:", value = 0,max=2,min=-2,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 150px;",
                               sliderInput("s_nlp_daa", "scale (mode = location ± scale * √2 ): ", value = .045,max=1,min=.045,ticks = FALSE,step = .001))
                         ),
                         
                       actionButton("run", label = "Run"),
                       p("Note: Error would occur if the required sample size is more than 100,000"),
                          
                       p("Please report any issues to the developer and the maintainer of the app, T.K. Wong, at the email address: t.k.wong3004@gmail.com ")),
                      
                       mainPanel(
                         fluidRow(
                           column(6, plotOutput("plot1")),
                           column(6,   htmlOutput("result"))
                         ),
                        plotOutput("plot2")
                       )
                       
                       )
              ),
              tabPanel("Independent t-test(equal variance)",
               sidebarLayout(sidebarPanel(
                 radioButtons("mode2", label = p(strong("BFDA:")), choices = c("for sample size determination","for a fixed sample size"), inline = TRUE),
                 p(strong("Specify Alternative hypothese")),
                 radioButtons("rb2", label = h4(HTML(paste0("H", tags$sub("1"), " :"))), choices = c("δ ≠ 0","δ > 0", "δ < 0"), inline = TRUE),
                

                 conditionalPanel(
                   condition = "input.rb2 == 'δ ≠ 0'||input.rb2 == 'δ > 0'||input.rb2 == 'δ < 0'",
                   radioButtons("model2", label = p(strong("Specify analysis prior for δ under H₁:")), choices = c("Cauchy","Normal", "t-student","Non-local"), inline = TRUE)),
                 
                 conditionalPanel(
                   condition = "(input.model2 == 'Cauchy') && (input.rb2 == 'δ ≠ 0'||input.rb2 == 'δ > 0'||input.rb2 == 'δ < 0')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_c2", "location:", value = 0,min=-2,max=2,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("s_c2", "scale:", value = .707,min=.1,max=3,ticks = FALSE,step = .001)),
                   
                 ),
                 conditionalPanel(
                   condition = "(input.model2 == 'Normal') && (input.rb2 == 'δ ≠ 0'||input.rb2 == 'δ > 0'||input.rb2 == 'δ < 0')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_n2", "Mean:", value = 0,max=2,min=-2,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("s_n2", "SD:", value = 1,max=3,min=.1,ticks = FALSE,step = .001)),
                   
                 ),
                 conditionalPanel(
                   condition = "(input.model2 == 't-student') && (input.rb2 == 'δ ≠ 0'||input.rb2 == 'δ > 0'||input.rb2 == 'δ < 0')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_t2", "location:", value = 0,max= 2,min=-2,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("s_t2", "scale: ", value = .707,min=.1,max=2,ticks = FALSE,step = .001)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("df_t2", "df: ", value = 1,min=1,max = 50,ticks = FALSE,step = 1))
                   
                 ),
                 conditionalPanel(
                   condition = "(input.model2 == 'Non-local') && (input.rb2 == 'δ ≠ 0'||input.rb2 == 'δ > 0'||input.rb2 == 'δ < 0')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_nlp2", "location:", value = 0,max=2,min=-2,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 150px;",
                       sliderInput("s_nlp2", "scale (mode = location ± scale * √2 ): ", value = .045,max=1,min=.045,ticks = FALSE,step = .001))
                 ),
                 conditionalPanel(
                   condition = "input.mode2 == 'for a fixed sample size'",
                   p(strong("Specify sample size per group")),
                   
                   div(style="display: inline-block; width: 100px;",
                       numericInput("n1", "Group 1:", value = 51, min = 2, max = 100000, step = 1)),
                   div(style="display: inline-block; width: 100px;",
                       numericInput("n2", "Group 2:", value = 50, min = 2, max = 100000, step = 1)),
                   div(style="display: inline-block; width: 100px;")

                   
                   ),
                 conditionalPanel(
                   condition = "input.mode2 == 'for sample size determination'",
                   sliderInput("r2", p(strong("Specify the ratio of sample sizes in two groups N2/N1:")), 
                                                                       min=1,max = 10,value = 1)
                   ),
                   

                 
                 sliderInput("de2", p(strong("Specify the bound of compelling evidence BF:")), min=1,max = 20,value = 3),
                 conditionalPanel(
                   condition = "input.mode2 == 'for sample size determination'",
                 sliderInput("pro2", p(strong("Specify the desired probability of true positive evidence")),
                             min =.5, max = .99, value = .8,step = .01)),
                 
                 radioButtons("daa2", label = p(strong("Design prior is the same as analysis prior :")), choices = c("Yes","No"), inline = TRUE),
                 
                 conditionalPanel(
                   condition = "(input.daa2 == 'No')",
                   radioButtons("model_daa2", label = p(strong("Specify the model for design prior:")), choices = c("Cauchy","Normal", "t-student","Non-local"), inline = TRUE)),
                 
                 conditionalPanel(
                   condition = "(input.model_daa2 == 'Cauchy') &&  (input.daa2 == 'No')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_c_daa2", "location:", value = 0,min=-2,max=2,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("s_c_daa2", "scale:", value = .707,min=.1,max=3,ticks = FALSE,step = .001)),
                   
                 ),
                 conditionalPanel(
                   condition = "(input.model_daa2 == 'Normal') && (input.daa2 == 'No')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_n_daa2", "Mean:", value = 0,max=2,min=-2,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("s_n_daa2", "SD:", value = 1,max=3,min=.1,ticks = FALSE,step = .001)),
                   
                 ),
                 conditionalPanel(
                   condition = "(input.model_daa2 == 't-student') && (input.daa2 == 'No')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_t_daa2", "location:", value = 0,max= 2,min=-2,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("s_t_daa2", "scale: ", value = .707,min=.1,max=2,ticks = FALSE,step = .001)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("df_t_daa2", "df: ", value = 1,min=1,max = 50,ticks = FALSE,step = 1))
                   
                 ),
                 conditionalPanel(
                   condition = "(input.model_daa2 == 'Non-local') && (input.daa2 == 'No')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_nlp_daa2", "location:", value = 0,max=2,min=-2,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 150px;",
                       sliderInput("s_nlp_daa2", "scale (mode = location ± scale * √2 ): ", value = .045,max=1,min=.045,ticks = FALSE,step = .001))
                 ),
          
                 
                 actionButton("run2", label = "Run"),
                 p("Note: Error would occur if the required sample size is more than 100,000"),
                 p("Please report any issues to the developer and the maintainer of the app, T.K. Wong, at the email address: t.k.wong3004@gmail.com ")
               ),
               mainPanel(
                 fluidRow(
                   column(6, plotOutput("plot12")),
                   column(6,  htmlOutput("result2"))
                 ),
                 plotOutput("plot22")
               )
               
               )
              )
              
              )
              





server <- function(input, output, session) {
  observeEvent(input$run,{
    probabilities_table <- reactive({
      hypothesis <-switch(input$rb,
                          "δ ≠ 0" = "!=",
                          "δ > 0" = ">",
                          "δ < 0" = "<"
      )
      model <- switch(input$model,
                      "Cauchy" = "Cauchy",
                      "Normal" = "Normal",
                      "t-student" = "t-distribution",
                      "Non-local" = "NLP"
      )
      location <- switch(model,
                         "Cauchy" = input$l_c,
                         "Normal" = input$l_n,
                         "t-distribution" = input$l_t,
                         "NLP" = input$l_nlp
      )
      scale <- switch(model,
                      "Cauchy" = input$s_c,
                      "Normal" = input$s_n,
                      "t-distribution" = input$s_t,
                      "NLP" = input$s_nlp
      )
      dff <- input$df_t
      
      N <- input$N
      mode <- switch(input$mode,
                     "for sample size determination" = 1,
                     "for a fixed sample size" = 0)
      hypothesis_d =hypothesis
      
      model_d <- switch(input$model_daa,
                      "Cauchy" = "Cauchy",
                      "Normal" = "Normal",
                      "t-student" = "t-distribution",
                      "Non-local" = "NLP"
      )
      location_d <- switch(model_d,
                         "Cauchy" = input$l_c_daa,
                         "Normal" = input$l_n_daa,
                         "t-distribution" = input$l_t_daa,
                         "NLP" = input$l_nlp_daa
      )
      scale_d <- switch(model_d,
                      "Cauchy" = input$s_c_daa,
                      "Normal" = input$s_n_daa,
                      "t-distribution" = input$s_t_daa,
                      "NLP" = input$s_nlp_daa
      )
      dff_d <- input$df_t_daa
      
      de_an_prior <- switch(input$daa,
                            'Yes' = 1,
                            'No' = 0
                            
      )
      
      
      D <- input$de
      target <-input$pro
      dat = suppressWarnings(Table(D,target,model,location ,scale,dff,hypothesis,model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior ,N,mode))
      if (location >2 | location< (-2)){
        dat = NA
      }
      dat
    })

    output$plot1 <- renderPlot({
      hypothesis <-switch(input$rb,
                          "δ ≠ 0" = "!=",
                          "δ > 0" = ">",
                          "δ < 0" = "<"
      )
      model <- switch(input$model,
                      "Cauchy" = "Cauchy",
                      "Normal" = "Normal",
                      "t-student" = "t-distribution",
                      "Non-local" = "NLP"
      )
      location <- switch(model,
                         "Cauchy" = input$l_c,
                         "Normal" = input$l_n,
                         "t-distribution" = input$l_t,
                         "NLP" = input$l_nlp
      )
      scale <- switch(model,
                      "Cauchy" = input$s_c,
                      "Normal" = input$s_n,
                      "t-distribution" = input$s_t,
                      "NLP" = input$s_nlp
      )
      dff <- input$df_t
      
      
      hypothesis_d =hypothesis
      model_d <- switch(input$model_daa,
                        "Cauchy" = "Cauchy",
                        "Normal" = "Normal",
                        "t-student" = "t-distribution",
                        "Non-local" = "NLP"
      )
      location_d <- switch(model_d,
                           "Cauchy" = input$l_c_daa,
                           "Normal" = input$l_n_daa,
                           "t-distribution" = input$l_t_daa,
                           "NLP" = input$l_nlp_daa
      )
      scale_d <- switch(model_d,
                        "Cauchy" = input$s_c_daa,
                        "Normal" = input$s_n_daa,
                        "t-distribution" = input$s_t_daa,
                        "NLP" = input$s_nlp_daa
      )
      dff_d <- input$df_t_daa
      
      de_an_prior <- switch(input$daa,
                            "Yes" = 1,
                            "No" = 0
                            
      )
      
      D <- input$de
      target <-input$pro
      suppressWarnings(prior_plot(D ,target ,model,location,scale,dff,hypothesis, model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior))
      
    })
    
    output$plot2 <- renderPlot({
      hypothesis <-switch(input$rb,
                          "δ ≠ 0" = "!=",
                          "δ > 0" = ">",
                          "δ < 0" = "<"
      )
      model <- switch(input$model,
                      "Cauchy" = "Cauchy",
                      "Normal" = "Normal",
                      "t-student" = "t-distribution",
                      "Non-local" = "NLP"
      )

      location <- switch(model,
                         "Cauchy" = input$l_c,
                         "Normal" = input$l_n,
                         "t-distribution" = input$l_t,
                         "NLP" = input$l_nlp,
                         "Point" = input$h1
      )
      scale <- switch(model,
                      "Cauchy" = input$s_c,
                      "Normal" = input$s_n,
                      "t-distribution" = input$s_t,
                      "NLP" = input$s_nlp,
                      "Point" = .707
      )
      dff <- input$df_t
      D <- input$de
      target <-input$pro
      df <-probabilities_table()
      if (location >2 | location < (-2)){
        plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10),axes = F,main = "The location parameter has a upper limit of 2 and a lower limit of -2")
      }else{
        suppressWarnings(bf10_t(D =D,df[1,6],target = target, model = model,location =location ,scale= scale,dff=dff,hypothesis =hypothesis))
        
      }
    })
    
    result <- reactive({
      dat <- probabilities_table()
      
      HTML(paste0(
        "<style>
      table {
        width: 80%;  /* Increase the width of the table */
        border-collapse: collapse;
        margin: 20px;
        font-family: 'Times New Roman', Times, serif;
      }
      th, td {
        border: 1px solid black;
        padding: 8px;
      }
      th {
        text-align: left;
        background-color: #f2f2f2;
      }
      .category th {
        border-bottom: 1px solid black; /* Add bottom border to category headers */
        padding-bottom: 10px; /* Add space below category headers */
      }
      .noborder {
        border: none !important; /* Remove all borders for cells with this class */
      }
    </style>",
        "<table>",
        "<tr><th colspan='2'>Probability of Compelling Evidence</th></tr>",
        "<tr><td class='noborder'>p(BF<sub>10</sub> > ", input$de, " | H<sub>1</sub>)</td>",
        "<td class='noborder'>", round(dat[1], 3), "</td></tr>",
        "<tr><td class='noborder'>p(BF<sub>01</sub> > ", input$de, " | H<sub>0</sub>)</td>",
        "<td class='noborder'>", round(dat[3], 3), "</td></tr>",
        "<tr><th colspan='2'>Probability of Misleading Evidence</th></tr>",
        "<tr><td class='noborder'>p(BF<sub>01</sub> > ", input$de, " | H<sub>1</sub>)</td>",
        "<td class='noborder'>", round(dat[2], 3), "</td></tr>",
        "<tr><td class='noborder'>p(BF<sub>10</sub> > ", input$de, " | H<sub>0</sub>)</td>",
        "<td class='noborder'>", round(dat[4], 3), "</td></tr>",
        "<tr><th colspan='2'>Required Sample Size</th></tr>",
        "<tr><td class='noborder'>N</td>",
        "<td class='noborder'>", dat[5], "</td></tr>",
        "<tr><td class='noborder'>Exact needed df</td>",
        "<td class='noborder'>", round(dat[6], 5), "</td></tr>",
        "</table>"
      ))
    })
    
    
    
    output$result <- renderUI({
      result()
    })
    
    
    
    
  })
  

  
  observeEvent(input$run2,{
    probabilities_table2 <- reactive({
      hypothesis <-switch(input$rb2,
                          "δ ≠ 0" = "!=",
                          "δ > 0" = ">",
                          "δ < 0" = "<"
      )
      model <- switch(input$model2,
                      "Cauchy" = "Cauchy",
                      "Normal" = "Normal",
                      "t-student" = "t-distribution",
                      "Non-local" = "NLP"
      )
      location <- switch(model,
                         "Cauchy" = input$l_c2,
                         "Normal" = input$l_n2,
                         "t-distribution" = input$l_t2,
                         "NLP" = input$l_nlp2
      )
      scale <- switch(model,
                      "Cauchy" = input$s_c2,
                      "Normal" = input$s_n2,
                      "t-distribution" = input$s_t2,
                      "NLP" = input$s_nlp2
      )
      dff <- input$df_t2
      hypothesis_d = hypothesis
      
      model_d <- switch(input$model_daa2,
                      "Cauchy" = "Cauchy",
                      "Normal" = "Normal",
                      "t-student" = "t-distribution",
                      "Non-local" = "NLP"
      )
      location_d <- switch(model_d,
                         "Cauchy" = input$l_c_daa2,
                         "Normal" = input$l_n_daa2,
                         "t-distribution" = input$l_t_daa2,
                         "NLP" = input$l_nlp_daa2
      )
      scale_d <- switch(model_d,
                      "Cauchy" = input$s_c_daa2,
                      "Normal" = input$s_n_daa2,
                      "t-distribution" = input$s_t_daa2,
                      "NLP" = input$s_nlp_daa2
      )
      dff_d <- input$df_t_daa2
      
      
      de_an_prior <- switch(input$daa2,
                            "Yes" = 1,
                            "No" = 0
      )
      mode2 <- switch(input$mode2,
                     "for sample size determination" = 1,
                     "for a fixed sample size" = 0)
      n1 = input$n1
      n2 = input$n2
      
      D <- input$de2
      target <-input$pro2
      r <- input$r2
      dat = suppressWarnings(Table_two(D,r,target,model,location,scale,dff, hypothesis, model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior,n1,n2,mode2))
      if (location >2 | location< (-2)){
        dat = NA
      }
      dat
    })
    
    output$plot12 <- renderPlot({
      hypothesis <-switch(input$rb2,
                          "δ ≠ 0" = "!=",
                          "δ > 0" = ">",
                          "δ < 0" = "<"
      )
      model <- switch(input$model2,
                      "Cauchy" = "Cauchy",
                      "Normal" = "Normal",
                      "t-student" = "t-distribution",
                      "Non-local" = "NLP"
      )

      location <- switch(model,
                         "Cauchy" = input$l_c2,
                         "Normal" = input$l_n2,
                         "t-distribution" = input$l_t2,
                         "NLP" = input$l_nlp2
      )
      scale <- switch(model,
                      "Cauchy" = input$s_c2,
                      "Normal" = input$s_n2,
                      "t-distribution" = input$s_t2,
                      "NLP" = input$s_nlp2
      )
      dff <- input$df_t2
      hypothesis_d = hypothesis
      model_d <- switch(input$model_daa2,
                        "Cauchy" = "Cauchy",
                        "Normal" = "Normal",
                        "t-student" = "t-distribution",
                        "Non-local" = "NLP"
      )
      location_d <- switch(model_d,
                           "Cauchy" = input$l_c_daa2,
                           "Normal" = input$l_n_daa2,
                           "t-distribution" = input$l_t_daa2,
                           "NLP" = input$l_nlp_daa2
      )
      scale_d <- switch(model_d,
                        "Cauchy" = input$s_c_daa2,
                        "Normal" = input$s_n_daa2,
                        "t-distribution" = input$s_t_daa2,
                        "NLP" = input$s_nlp_daa2
      )
      dff_d <- input$df_t_daa2
      

      
      de_an_prior <- switch(input$daa2,
                            "Yes" = 1,
                            "No" = 0
                            
      )
      
      
      
      
      D <- input$de2
      target <-input$pro2
      
      
      suppressWarnings(prior_plot(D ,target ,model,location,scale,dff,hypothesis,model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior ))
      
      
    })
    
    output$plot22 <- renderPlot({
      hypothesis <-switch(input$rb2,
                          "δ ≠ 0" = "!=",
                          "δ > 0" = ">",
                          "δ < 0" = "<"
      )
      model <- switch(input$model2,
                      "Cauchy" = "Cauchy",
                      "Normal" = "Normal",
                      "t-student" = "t-distribution",
                      "Non-local" = "NLP"
      )
      
      location <- switch(model,
                         "Cauchy" = input$l_c2,
                         "Normal" = input$l_n2,
                         "t-distribution" = input$l_t2,
                         "NLP" = input$l_nlp2
      )
      scale <- switch(model,
                      "Cauchy" = input$s_c2,
                      "Normal" = input$s_n2,
                      "t-distribution" = input$s_t2,
                      "NLP" = input$s_nlp2
      )
      dff <- input$df_t2
      D <- input$de2
      target <-input$pro2
      ddf <-probabilities_table2()
      
      if (location >2 | location < (-2)){
        plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10),axes = F,main = "The location parameter has a upper limit of 2 and a lower limit of -2")
      }else{
        suppressWarnings(bf10_two(D ,ddf[1,5],input$r2, target,model,location ,scale,dff, hypothesis ))
        
      }
    })
    
    
   result32 <- reactive({
  dat <- probabilities_table2()
  
  HTML(paste0(
    "<style>
      table {
        width: 80%;  /* Increase the width of the table */
        border-collapse: collapse;
        margin: 20px;
        font-family: 'Times New Roman', Times, serif;
      }
      th, td {
        border: 1px solid black;
        padding: 8px;
      }
      th {
        text-align: left;
        background-color: #f2f2f2;
      }
      .category th {
        border-bottom: 1px solid black; /* Add bottom border to category headers */
        padding-bottom: 10px; /* Add space below category headers */
      }
      .noborder {
        border: none !important; /* Remove all borders for cells with this class */
      }
    </style>",
    "<table>",
    "<tr><th colspan='2'>Probability of Compelling Evidence</th></tr>",
    "<tr><td class='noborder'>p(BF<sub>10</sub> > ", input$de2, " | H<sub>1</sub>)</td>",
    "<td class='noborder'>", round(dat[1], 3), "</td></tr>",
    "<tr><td class='noborder'>p(BF<sub>01</sub> > ", input$de2, " | H<sub>0</sub>)</td>",
    "<td class='noborder'>", round(dat[3], 3), "</td></tr>",
    "<tr><th colspan='2'>Probability of Misleading Evidence</th></tr>",
    "<tr><td class='noborder'>p(BF<sub>01</sub> > ", input$de2, " | H<sub>1</sub>)</td>",
    "<td class='noborder'>", round(dat[2], 3), "</td></tr>",
    "<tr><td class='noborder'>p(BF<sub>10</sub> > ", input$de2, " | H<sub>0</sub>)</td>",
    "<td class='noborder'>", round(dat[4], 3), "</td></tr>",
    "<tr><th colspan='2'>Required Sample Size</th></tr>",
    "<tr><td class='noborder'>Group 1</td>",
    "<td class='noborder'>", dat[5], "</td></tr>",
    "<tr><td class='noborder'>Group 2</td>",
    "<td class='noborder'>", round(ceiling(dat[6]), 0), "</td></tr>",
    "</table>"
  ))
})
 
    output$result2 <- renderUI({
      result32()
    })
    

  })
  }
shinyApp(ui, server)


