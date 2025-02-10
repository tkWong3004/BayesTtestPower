library(rootSolve)
library(shiny)
library(Rcpp)
library(BH)
library("bslib")
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
                                 sliderInput("l_c", "location:", value = 0,min=-3,max=3,ticks = FALSE,step = .01)),
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
                               sliderInput("l_t", "location:", value = 0,max= 3,min=-3,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("s_t", "scale: ", value = .707,min=.1,max=2,ticks = FALSE,step = .001)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("df_t", "df: ", value = 1,min=1,max = 50,ticks = FALSE,step = 1))
                           
                         ),
                         conditionalPanel(
                           condition = "(input.model == 'Non-local') && (input.rb == 'δ ≠ 0'||input.rb == 'δ > 0'||input.rb == 'δ < 0')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_nlp", "location:", value = 0,max=3,min=-3,ticks = FALSE,step = .01)),
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
                           radioButtons("model_daa", label = p(strong("Specify the model for design prior:")), choices = c("Cauchy","Normal", "t-student","Non-local","Point"), inline = TRUE)),
                         
                         conditionalPanel(
                           condition = "(input.model_daa == 'Cauchy') &&  (input.daa == 'No')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_c_daa", "location:", value = 0,min=-3,max=3,ticks = FALSE,step = .01)),
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
                               sliderInput("l_t_daa", "location:", value = 0,max= 3,min=-3,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("s_t_daa", "scale: ", value = .707,min=.1,max=2,ticks = FALSE,step = .001)),
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("df_t_daa", "df: ", value = 1,min=1,max = 50,ticks = FALSE,step = 1))
                           
                         ),
                         conditionalPanel(
                           condition = "(input.model_daa == 'Non-local') && (input.daa == 'No')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("l_nlp_daa", "location:", value = 0,max=3,min=-3,ticks = FALSE,step = .01)),
                           div(style="display: inline-block; width: 150px;",
                               sliderInput("s_nlp_daa", "scale (mode = location ± scale * √2 ): ", value = .045,max=1,min=.045,ticks = FALSE,step = .001))
                         ),
                         conditionalPanel(
                           condition = "(input.model_daa == 'Point') && (input.daa == 'No') && (input.rb == 'δ ≠ 0')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("p2t", "location:", value = 0,max=3,min=-3,ticks = FALSE,step = .01)),
                         ),conditionalPanel(
                           condition = "(input.model_daa == 'Point') && (input.daa == 'No') && (input.rb == 'δ > 0')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("p1tg", "location:", value = .01,max=3,min=.01,ticks = FALSE,step = .01)),
                         ),
                         conditionalPanel(
                           condition = "(input.model_daa == 'Point') && (input.daa == 'No') && (input.rb == 'δ < 0')",
                           div(style="display: inline-block; width: 100px;",
                               sliderInput("p1tl", "location:", value = -.01,max=-.01,min=-3,ticks = FALSE,step = .01)),
                         ),
                       actionButton("run", label = "Run"),
                       p("Note: Error would occur if the required sample size is more than 100,000"),
                          
                       p("Please report any issues to the developer and the maintainer of the app, T.K. Wong, at the email address: t.k.wong3004@gmail.com ")),

                       mainPanel(
                         navset_card_underline(
                           nav_panel("Result",
                         fluidRow(
                           column(6, plotOutput("plot1")),
                           column(6,   htmlOutput("result"))
                         ),
                         plotOutput("plot2")),
                         nav_panel("Power Curve",plotOutput("plot3"))
                         
                         )
                         
                         
                       )
                         
                       )
                       
                       
                       
              ),
              tabPanel("Independent samples t-test(equal variance)",
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
                       sliderInput("l_c2", "location:", value = 0,min=-3,max=3,ticks = FALSE,step = .01)),
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
                       sliderInput("l_t2", "location:", value = 0,max= 3,min=-3,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("s_t2", "scale: ", value = .707,min=.1,max=2,ticks = FALSE,step = .001)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("df_t2", "df: ", value = 1,min=1,max = 50,ticks = FALSE,step = 1))
                   
                 ),
                 conditionalPanel(
                   condition = "(input.model2 == 'Non-local') && (input.rb2 == 'δ ≠ 0'||input.rb2 == 'δ > 0'||input.rb2 == 'δ < 0')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_nlp2", "location:", value = 0,max=3,min=-3,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 150px;",
                       sliderInput("s_nlp2", "scale (mode = location ± scale * √2 ): ", value = .045,max=1,min=.045,ticks = FALSE,step = .001))
                 ),
                 conditionalPanel(
                   condition = "input.mode2 == 'for a fixed sample size'",
                   p(strong("Specify sample size per group")),
                   
                   div(style="display: inline-block; width: 100px;",
                       numericInput("n1", "Group 1:", value = 52, min = 2, max = 100000, step = 1)),
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
                   radioButtons("model_daa2", label = p(strong("Specify the model for design prior:")), choices = c("Cauchy","Normal", "t-student","Non-local","Point"), inline = TRUE)),
                 
                 conditionalPanel(
                   condition = "(input.model_daa2 == 'Cauchy') &&  (input.daa2 == 'No')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_c_daa2", "location:", value = 0,min=-3,max=3,ticks = FALSE,step = .01)),
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
                       sliderInput("l_t_daa2", "location:", value = 0,max= 3,min=-3,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("s_t_daa2", "scale: ", value = .707,min=.1,max=2,ticks = FALSE,step = .001)),
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("df_t_daa2", "df: ", value = 1,min=1,max = 50,ticks = FALSE,step = 1))
                   
                 ),
                 conditionalPanel(
                   condition = "(input.model_daa2 == 'Non-local') && (input.daa2 == 'No')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("l_nlp_daa2", "location:", value = 0,max=3,min=-3,ticks = FALSE,step = .01)),
                   div(style="display: inline-block; width: 150px;",
                       sliderInput("s_nlp_daa2", "scale (mode = location ± scale * √2 ): ", value = .045,max=1,min=.045,ticks = FALSE,step = .001))
                 ),
                 conditionalPanel(
                   condition = "(input.model_daa2 == 'Point') && (input.daa2 == 'No') && (input.rb2 == 'δ ≠ 0')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("p2t2", "location:", value = 0,max=3,min=-3,ticks = FALSE,step = .01)),
                 ),conditionalPanel(
                   condition = "(input.model_daa2 == 'Point') && (input.daa2 == 'No') && (input.rb2 == 'δ > 0')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("p1tg2", "location:", value = .01,max=3,min=.01,ticks = FALSE,step = .01)),
                 ),
                 conditionalPanel(
                   condition = "(input.model_daa2 == 'Point') && (input.daa2 == 'No') && (input.rb2 == 'δ < 0')",
                   div(style="display: inline-block; width: 100px;",
                       sliderInput("p1tl2", "location:", value = -.01,max=-.01,min=-3,ticks = FALSE,step = .01)),
                 ),
          
                 
                 actionButton("run2", label = "Run"),
                 p("Note: Error would occur if the required sample size is more than 100,000"),
                 p("Please report any issues to the developer and the maintainer of the app, T.K. Wong, at the email address: t.k.wong3004@gmail.com ")
               ),
               
               mainPanel(
                 
                 
                 navset_card_underline(
                   nav_panel("Result",
                 fluidRow(
                   column(6, plotOutput("plot12")),
                   column(6,  htmlOutput("result2"))
                 ),plotOutput("plot22")),
                 nav_panel("Power Curve",plotOutput("plot32")))
                 
                 
                 
                 )
               
               )
              )
              
              )
              





server <- function(input, output, session) {
  input_t1 <- reactive({
    mode <- switch(input$mode,
                   "for sample size determination" = 1,
                   "for a fixed sample size" = 0)
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

    hypothesis_d =hypothesis
    
    model_d <- switch(input$model_daa,
                      "Cauchy" = "Cauchy",
                      "Normal" = "Normal",
                      "t-student" = "t-distribution",
                      "Non-local" = "NLP",
                      "Point" = "Point"
    )
    if (model_d != "Point"){
      location_d <- switch(model_d,
                           "Cauchy" = input$l_c_daa,
                           "Normal" = input$l_n_daa,
                           "t-distribution" = input$l_t_daa,
                           "NLP" = input$l_nlp_daa
      )}else{
        location_d <-switch(hypothesis,
                            "!=" = input$p2t,
                            ">" = input$p1tg,
                            "<" = input$p1tl)
      }
    scale_d <- switch(model_d,
                      "Cauchy" = input$s_c_daa,
                      "Normal" = input$s_n_daa,
                      "t-distribution" = input$s_t_daa,
                      "NLP" = input$s_nlp_daa,
                      "Point" = 0
    )
    dff_d <- input$df_t_daa
    
    de_an_prior <- switch(input$daa,
                          'Yes' = 1,
                          'No' = 0
                          
    )
    
    
    DD <- input$de
    target <-input$pro 
    
    list(
      hypothesis = hypothesis,
      model = model,
      location = location,
      scale = scale,
      dff = dff,
      N = N,
      mode = mode,
      hypothesis_d =hypothesis,
      model_d = model_d,
      location_d = location_d,
      scale_d = scale_d,
      dff_d = dff_d,
      de_an_prior = de_an_prior,
      DD = DD,
      target = target
    )
    
    
  })
  
  
  observeEvent(input$run,{
    xx = input_t1()
    
    dat = suppressWarnings(Table(xx$DD,
                                   xx$target,
                                   xx$model,
                                   xx$location ,
                                   xx$scale,
                                   xx$dff,
                                   xx$hypothesis,
                                   xx$model_d,
                                   xx$location_d,
                                   xx$scale_d,
                                   xx$dff_d, 
                                   xx$hypothesis_d,
                                   xx$de_an_prior ,
                                   xx$N,
                                   xx$mode))

    output$plot1 <- renderPlot({
      suppressWarnings(prior_plot(xx$DD ,
                                  xx$target ,
                                  xx$model,
                                  xx$location,
                                  xx$scale,
                                  xx$dff,
                                  xx$hypothesis, 
                                  xx$model_d,
                                  xx$location_d,
                                  xx$scale_d,
                                  xx$dff_d, 
                                  xx$hypothesis_d,
                                  xx$de_an_prior))
      
    })
    
    output$plot2 <- renderPlot({
        suppressWarnings(bf10_t(D =xx$DD,dat[1,6],target = xx$target, model = xx$model,
                                location =xx$location ,scale= xx$scale,dff=xx$dff,hypothesis =xx$hypothesis))
        
    })

    
    output$result <- renderUI({
  
      
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
        "<tr><td class='noborder'>p(BF<sub>10</sub> > ", xx$DD, " | H<sub>1</sub>)</td>",
        "<td class='noborder'>", round(dat[1], 3), "</td></tr>",
        "<tr><td class='noborder'>p(BF<sub>01</sub> > ", xx$DD, " | H<sub>0</sub>)</td>",
        "<td class='noborder'>", round(dat[3], 3), "</td></tr>",
        "<tr><th colspan='2'>Probability of Misleading Evidence</th></tr>",
        "<tr><td class='noborder'>p(BF<sub>01</sub> > ", xx$DD, " | H<sub>1</sub>)</td>",
        "<td class='noborder'>", round(dat[2], 3), "</td></tr>",
        "<tr><td class='noborder'>p(BF<sub>10</sub> > ", xx$DD, " | H<sub>0</sub>)</td>",
        "<td class='noborder'>", round(dat[4], 3), "</td></tr>",
        "<tr><th colspan='2'>Required Sample Size</th></tr>",
        "<tr><td class='noborder'>N</td>",
        "<td class='noborder'>", dat[5], "</td></tr>",
        "<tr><td class='noborder'>Exact needed df</td>",
        "<td class='noborder'>", round(dat[6], 5), "</td></tr>",
        "</table>"
      ))
    })
    

    output$plot3 <- renderPlot({
      Power_t1(xx$DD,xx$model,xx$location,xx$scale,xx$dff, xx$hypothesis, xx$model_d,xx$location_d,xx$scale_d,xx$dff_d, xx$de_an_prior,dat[1,5],xx$target ,xx$mode) 
      })

    
  })
  
  input_t2 <- reactive({
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
                      "Non-local" = "NLP",
                      "Point" = "Point"
    )
    if (model_d != "Point"){
      location_d <- switch(model_d,
                           "Cauchy" = input$l_c_daa2,
                           "Normal" = input$l_n_daa2,
                           "t-distribution" = input$l_t_daa2,
                           "NLP" = input$l_nlp_daa2
      )}else{
        location_d <-switch(hypothesis,
                            "!=" = input$p2t2,
                            ">" = input$p1tg2,
                            "<" = input$p1tl2)
      }
    
    
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
    
    list(
      hypothesis = hypothesis,
      model = model,
      location = location,
      scale = scale,
      dff = dff,
      hypothesis_d = hypothesis,
      mode2  =mode2 ,
      hypothesis_d =hypothesis,
      model_d = model_d,
      location_d = location_d,
      scale_d = scale_d,
      dff_d = dff_d,
      de_an_prior = de_an_prior,
      D = D,
      target = target,
      n1= n1,
      n2=n2,
      r=r
    )
    
    
  })
  
  observeEvent(input$run2,{
    x=input_t2()

      
      dat2 <- suppressWarnings(Table_two(x$D,x$r,x$target,x$model,x$location,x$scale,
                                       x$dff, x$hypothesis, 
                                       x$model_d,x$location_d,
                                       x$scale_d,x$dff_d, x$hypothesis_d,
                                       x$de_an_prior,x$n1,x$n2,x$mode2))
  
    
    output$plot12 <- renderPlot({
      
      suppressWarnings(prior_plot(x$D ,x$target ,x$model,x$location,x$scale,
                                  x$dff,x$hypothesis,x$model_d,x$location_d,
                                  x$scale_d,x$dff_d, x$hypothesis_d,x$de_an_prior ))
      
      
    })
    
    output$plot22 <- renderPlot({

        suppressWarnings(bf10_two(x$D ,dat2[1,5],x$r, x$target,x$model,
                                  x$location ,x$scale,x$dff, x$hypothesis ))
        
      
    })
    
    
   result32 <- reactive({
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
    "<tr><td class='noborder'>p(BF<sub>10</sub> > ", x$D, " | H<sub>1</sub>)</td>",
    "<td class='noborder'>", round(dat2[1], 3), "</td></tr>",
    "<tr><td class='noborder'>p(BF<sub>01</sub> > ", x$D, " | H<sub>0</sub>)</td>",
    "<td class='noborder'>", round(dat2[3], 3), "</td></tr>",
    "<tr><th colspan='2'>Probability of Misleading Evidence</th></tr>",
    "<tr><td class='noborder'>p(BF<sub>01</sub> > ", x$D, " | H<sub>1</sub>)</td>",
    "<td class='noborder'>", round(dat2[2], 3), "</td></tr>",
    "<tr><td class='noborder'>p(BF<sub>10</sub> > ", x$D, " | H<sub>0</sub>)</td>",
    "<td class='noborder'>", round(dat2[4], 3), "</td></tr>",
    "<tr><th colspan='2'>Required Sample Size</th></tr>",
    "<tr><td class='noborder'>Group 1</td>",
    "<td class='noborder'>", dat2[5], "</td></tr>",
    "<tr><td class='noborder'>Group 2</td>",
    "<td class='noborder'>", round(ceiling(dat2[6]), 0), "</td></tr>",
    "</table>"
  ))
})
 
    output$result2 <- renderUI({
      result32()
    })
    
    output$plot32 <- renderPlot({
      Power_t2(x$D,x$r,x$model,x$location,x$scale,x$dff, x$hypothesis, x$model_d,x$location_d,x$scale_d,x$dff_d, x$de_an_prior,dat2[1,5],dat2[1,6],x$mode2 ,x$target) 
    })
   
  })
  }
shinyApp(ui, server)


