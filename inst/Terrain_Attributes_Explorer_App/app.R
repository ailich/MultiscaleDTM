kns2<- function(a,b,c,d,e){
    #Profc: Normal Slope Line Curvature
    out<- (-2 * (a*d^2 + c*d*e + b*e^2)) / ((d^2 + e^2)*(1 + d^2 + e^2)^1.5)
    return(out)
}

knc2<- function(a,b,c,d,e){
    #Planc: Normal Contour Curvature
    out<- -2*(a*(e^2) - c*d*e + b*d^2)/((d^2+e^2) * sqrt(1+d^2+e^2))
    return(out)
}

tgc2<- function(a,b,c,d,e){
    #TwistC: Contour geodesic torsion
    out<- (2*d*e*(a-b) - c*(d^2-e^2))/((d^2+e^2)*(1+d^2+e^2))
    return(out)
}

kmean2<- function(a,b,c,d,e){
    #Mean Curvature
    out<- -(a*(1+e^2) - c*d*e + b*(1+d^2)) / (sqrt((1+d^2+e^2)^3))
    return(out)
}

ku2<- function(a,b,c,d,e){
    #unsphericity curvature
    out<- sqrt(((a*(1+e^2) - c*d*e +b*(1+d^2)) / (sqrt((1+d^2+e^2)^3)))^2 - ((4*a*b-c^2)/(1+d^2+e^2)^2))
    return(out)
}

kmin2<- function(a,b,c,d,e){
    #Min Curvature
    out<- kmean2(a,b,c,d,e)-ku2(a,b,c,d,e)
    return(out)
}

kmax2<- function(a,b,c,d,e){
    #Max Curvature
    out<- kmean2(a,b,c,d,e)+ku2(a,b,c,d,e)
    return(out)
}

classify_features2<- function(slope, planc, maxc, minc) {
    slope_tolerance<- 1
    curvature_tolerance<- 0.0001
    dplyr::case_when(is.na(slope) ~ NA_character_,
                 (slope > slope_tolerance) & (planc > curvature_tolerance) ~ "Ridge",
                 (slope > slope_tolerance) & (planc < -curvature_tolerance) ~ "Channel",
                 slope > slope_tolerance ~ "Planar Slope",
                 (maxc > curvature_tolerance) & (minc > curvature_tolerance) ~ "Peak",
                 (maxc > curvature_tolerance) & (minc < -curvature_tolerance) ~ "Pass",
                 maxc > curvature_tolerance ~ "Ridge",
                 (minc < -curvature_tolerance) & (maxc < -curvature_tolerance) ~ "Pit",
                 minc < -curvature_tolerance ~ "Channel",
                 TRUE ~ "Planar Flat")
    }

nr<-3
nc<- 3
x<- matrix(1:nc, nrow=nr, ncol=nc, byrow=TRUE)
x<- x-mean(x)

y<- matrix(nr:1, nrow=nr, ncol=nc, byrow=FALSE)
y<- y-mean(y)

plot_surface<- function(a,b,c,d,e){
    z<- a*x^2+ b*y^2+ c*x*y + d*x + e*y
    df<- data.frame(x=as.vector(x),y=as.vector(y),z=as.vector(z))
    m<- lm(z ~ I(x^2)+I(y^2)+I(x*y)+x+y, data = df)
    if(any(round(m$coefficients, 1)[c(2,3,4,5,6,1)]!=c(a,b,c,d,e,0))){
        stop("Error: fit regression parameters don't match specified parameters")
    }
    
    asp<- (-pi/2) - atan2(e,d) #Shift aspect so north is zero
    asp[asp < 0]<- asp[asp < 0] + 2*pi
    asp[asp >= 2*pi]<- asp[asp >= 2*pi] - 2*pi # Constrain aspect from 0 to 2pi
    eastness<- sin(asp)
    northness<- cos(asp)

    rgl::plot3d(m)
    if(any(e!=0, d!=0)){
        rgl::arrow3d(p0=c(0,0,0), p1 = c(eastness, northness,0))
    }
}

make_table<- function(a,b,c,d,e){
    slp<- atan(sqrt(d^2 + e^2))
    asp<- (-pi/2) - atan2(e,d) #Shift aspect so north is zero
    asp[asp < 0]<- asp[asp < 0] + 2*pi
    asp[asp >= 2*pi]<- asp[asp >= 2*pi] - 2*pi # Constrain aspect from 0 to 2pi
    if(slp==0){
        asp<- NA_real_
        }
    eastness<- sin(asp)
    northness<- cos(asp)
    asp_deg<- asp * (180/pi)
    slp_deg<- slp * (180/pi)
    
    profc<- kns2(a,b,c,d,e)
    planc<- knc2(a,b,c,d,e)
    twistc<- tgc2(a,b,c,d,e)
    maxc<- kmax2(a,b,c,d,e)
    minc<- kmin2(a,b,c,d,e)
    meanc<- kmean2(a,b,c,d,e)
    
    #Translation to Minar 2020 formulas
    # zx<- (2*a*x) + (c*y) + d
    # zxx<- 2*a
    # zy<- (2*b*y) + (c*x) + e
    # zyy<- 2*b
    # zxy<- c
    
    #At central cell x and y are 0. Therefore,
    # zx<- d
    # zxx<- 2*a
    # zy<- e
    # zyy<- 2*b
    # zxy<- c
    
    out<- round(data.frame(profc = profc,
                           planc= planc,
                           twistc=twistc,
                           meanc=meanc,
                           maxc=maxc,
                           minc=minc,
                           slope_degrees=slp_deg,
                           aspect_degrees=asp_deg,
                           eastness=eastness,
                           northness=northness,
                           feature = NA_real_),3)
    
    out$feature[1] = classify_features2(slp_deg, planc, maxc, minc)
    return(out)
    }

ui <- fluidPage(

    # Application title
    titlePanel("Terrain Attributes Explorer"),

    # Sidebar with a slider input for number of bins 
    fluidRow(column(2,
        wellPanel("Regression Coefficients",
            sliderInput(inputId = "a", label = "a", min = -3, max = 3, value = 0, step = 0.1),
            sliderInput(inputId = "b", label = "b", min = -3, max = 3, value = 0, step = 0.1),
            sliderInput(inputId = "c", label = "c", min = -3, max = 3, value = 0, step = 0.1),
            sliderInput(inputId = "d", label = "d", min = -3, max = 3, value = 0, step = 0.1),
            sliderInput(inputId = "e", label = "e", min = -3, max = 3, value = 0, step = 0.1))), 
        column(5,rgl::rglwidgetOutput("plot",  width = 800, height = 600)),
        column(2, wellPanel(p("The shape of a small area of a digital elevation model can be approximated with a quadratic formulated as:"),
                            p(withMathJax("$$\\small{Z=aX^2+bY^2+cXY+dX+eY+f}$$")),
                            p("where Z is the elevation or depth, X is east/west direction, Y is north/south direction, and a-e are regression coefficients describing the shape of the quadratic surface and f is the intercept. From this surface we can calculate various terrain attributes including slope, aspect, and measures of curvature. Aspect is the compass direction of the slope measured as degrees clockwise from North, and can be decomposed into its north/south (northness) and east/west (eastness) components. There are three basic types of curvature: profile curvature which measures curvature along the direction of maximum slope, plan curvature which measures curvature perpendicular to the direction of maximum slope, and twisting curvature."))),
        column(3, wellPanel(p(withMathJax("$$\\small{\\text{Profile Curvature} = \\frac{-2(ad^2 + cde + be^2)}{(d^2+e^2)(1 + d^2 + e^2)^\\frac{3}{2}}}$$"),
                            p(withMathJax("$$\\small{\\text{Plan Curvature} = \\frac{-2(ae^2 - cde + bd^2)}{(d^2+e^2) \\sqrt{1+d^2+e^2}}}$$")),
                            p(withMathJax("$$\\small{\\text{Twisting Curvature} = \\frac{2de(a-b) - c(d^2-e^2)}{(d^2+e^2)(1+d^2+e^2)}}$$")),
                            p(withMathJax("$$\\small{\\text{Mean Curvature} = -\\frac{a(1+e^2) - cde + b(1+d^2)}{(1+d^2+e^2)^\\frac{3}{2}}}$$")),
                            p(withMathJax("$$\\small{\\text{ku} =\\sqrt{\\left(\\frac{a(1+e^2) - cde + b(1+d^2)}{(1+d^2+e^2)^\\frac{3}{2}}\\right)^2 - \\frac{4ab-c^2}{(1+d^2+e^2)^2}}}$$")),
                            p(withMathJax("$$\\small{\\text{Max Curvature} = \\text{Mean Curvature} + \\text{ku}}$$")),
                            p(withMathJax("$$\\small{\\text{Min Curvature} = \\text{Mean Curvature} - \\text{ku}}$$")),
                            p(withMathJax("$$\\small{\\text{Slope} = arctan(\\sqrt{d^2 + e^2})}$$")),
                            p(withMathJax("$$\\small{\\text{Aspect} = -\\frac{\\pi}{2} - arctan2(e,d)}$$")),
                            p(withMathJax("$$\\small{\\text{Eastness} = sin(\\text{Aspect})}$$")),
                            p(withMathJax("$$\\small{\\text{Northness} = cos(\\text{Aspect})}$$")),
                            ))),
        fluidRow(column(12, wellPanel(tableOutput("table"))))))

# Define server logic required to draw a histogram

server <- function(input, output, output2) {
    output$plot <- rgl::renderRglwidget({
        plot_surface(input$a, input$b, input$c, input$d, input$e)
        rgl::rglwidget()
    })

    output$table<- renderTable({make_table(input$a, input$b, input$c, input$d, input$e)}, bordered=TRUE)
    }

# Run the application 
shinyApp(ui = ui, server = server)
