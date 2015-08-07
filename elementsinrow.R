library(shiny)

# Complement generator function ####

compgen=function(oligo1,...) {
  
  # This function returns the PURE DNA complemnt of any single strand (sequence) handled as a string "atgctgctgcatc" written 5'-3' 
  #OLIGO1 MUST BE WRITTEN 5'-3'  
  
  #error checking example  
  #oligo1="atgcgcat"
  #print(length(oligo1)) Check
  #print(oligo1)
  
  # remove LNAS from sequence RETURNS ONLY THE PURE DNA COMPLEMENT
  oligo1=gsub(pattern="\\+","",oligo1)
  oligo1=tolower(oligo1)
  #print(oligo1) #check
  
  
  #expand the sequence
  sequence=substring(oligo1,1:nchar(oligo1), 1:nchar(oligo1)) # turns the string into a charcter vector eg. [1] "a" "t" g" "g" "g" "c"
  #print(sequence) CHECK
  
  
  
  
  
  complement=c(rep(0,length(sequence))) # make empty vector the same length as complement
  
  
  ## FILL THE VECTOR IN BACKWARDS SO WE DONT NEED TO FLIP IT AFTER
  #print(length(sequence)) CHECK
  
  for (i in 1:length(sequence)){ 
    
    #print(i)# CHECK
    base=sequence[i]
    
    if (identical(base,"a")){
      
      complement[length(complement)-i+1]="t"
      
    }
    else if (identical (base,"t")){
      complement[length(complement)-i+1]="a"
    }
    
    else if (identical(base,"g")){
      complement[length(complement)-i+1]="c"
    }
    else if (identical(base,"c")){
      complement[length(complement)-i+1]="g"
    }
    
    
    
  } # end of for loop
  
  
  ## now need to reverse the orders of the letters 
  
  
  
  ##  collapse into a single string
  
  #print(complement)# CHECK
  complement=paste(complement,sep="", collapse="")
  #print(complement)# CHECK
  #print(oligo1)
  
  
  return(complement)
  
}

# LIST GENERATOR FUNCTION ####

getComplementList=function (base){
  
  #return a list of choices that include all except the complementary base
  #input can be +A, +a, a, or A
  #base="c"
  
  
  #need to condition LNA bases to get the actual base
  if (nchar(base)>1){
    base=substring(base,2)
  }
  
  #turn to upper case
  base=toupper(base)
  
  #generate list
  
  if (base=="A"){
    compList=list("A","T","G","C")
  }
  else if (base=="T"){
    compList=list("T","A","G","C")
  }
  else if (base=="G"){
    compList=list("G","A","T","C")
  }
  else if (base=="C"){
    compList=list("C","A","T","G")
  }
  
  
  return(compList)
  
}

# ISOLATE LNA ####

isolate_LNA =function(oligo){
  
  #expand oligo
  
  oligo1=substring(oligo,1:nchar(oligo),1:nchar(oligo)) # expand it into a vector
  
  
  # inital loop params
  i=1
  oligo1mod=NULL
  while (i <= length(oligo1)){
    
    
    if(oligo1[i]=="+"){
      plusLNA=c(oligo1[i],oligo1[i+1])
      oligo1mod=c(oligo1mod,paste(plusLNA,sep="",collapse=""))
      
      i=i+2}
    
    
    
    else {
      oligo1mod=c(oligo1mod,oligo1[i])
      i=i+1
    }
    
  } #CONCATENATE "+" , "c" to "+c" when applicable so we can scan better
  
  
  return(oligo1mod)
  
  
}

# SEQ CHECK FUNCTION (INPUT CONDITIONING) ####

seqcheck=function(oligo1){
  
  # this function is to ensure that only the proper inputs get passed on to the tm pred program.
  # IT ALSO CONDITIONS THE OLIGO TO LOWER CASE 
  
  
  #oligo1= "+a+a+gc" Example oligo
  
  ### Preliminary Conditioning ###
  
  #remove spaces
  oligo1=gsub(" ","",oligo1)
  
  #put it to lower case
  oligo1=tolower(oligo1)
  
  
  
  
  # ANY OF THE FOLLOWING SEQUENCE MISTAKES WILL GIVE AN ERROR #
  nchar(oligo1)
  
  #0. Make sure there is an input
  if (nchar(oligo1)==0 |is.na(oligo1)){
    stop("Please enter a sequence written 5' to 3'")
  }
  #print("passed check 0") 
  

  
  #1.NON Watson CRICK or + checaters (we have some characters)
  
  nonwc=grep(pattern="[^atgc+ATGC]",oligo1,value=TRUE)
  
  if(length(nonwc)>0){ # the oligo contains non W/C or + characters
    
    stop("Only a,t,g,c bases and + are allowed")
    
  }
  
  
  #print("passed check 1")

  #2. CHECK THAT + is always followed by a letter (a,t,g,c)
  
  oligo1exp=substring(oligo1,1:nchar(oligo1),1:nchar(oligo1))
  for (i in 1:(length(oligo1exp))){
    
    if(oligo1exp[i]=="+"){# found a +
      # if the next element isn't a letter, give error
      #print(oligo1exp[i])
      #print(oligo1exp[i+1])
      
      nextletter=grep(pattern="[^atgcATGC]",oligo1exp[i+1], value=TRUE) # vector with non atgc elements
      
      if(length(nextletter)>0){
        stop("a + must be followed by an a,t,g, or c base")
        
      }
      
      if(is.na(oligo1exp[i+1])){
        stop("a + must be followed by an a,t,g, or c base")
        
      }
      
      
      
      
    }
  } # end of for loop
  
  #print("passed check 2")
  

  
  #3. Check that the duplex isnt too short (4 bp is shortest allowed _ AS PER IDT BIOPHYSICS)
  
  oligo1short=gsub("\\+", "", oligo1) # get rid of all +s 
  
  if(nchar(oligo1short)<4 & nchar(oligo1short)>0){
    
    stop ("sequence is too short")
    
  }
  
  #print("passed check 3")
  
  return(oligo1)
  
  
}

#REVERSE CHARACTERS ####

reverse_chars=function(string){
  
  string_split=strsplit(as.character(string),split="")
  reversed_split=string_split[[1]][nchar(string):1]
  paste(reversed_split,collapse="")
  
}

#  CSS   CHANGES  ####

css <- "
#sequence {
width: 500px;
}
#mismatches {
white-space: nowrap;
overflow-x: auto;
}
#mismatches .sequence-location {
white-space: normal;
display: inline-block;
width: 50px;
margin-right: 5px;
text-align: center;
}
#mismatches .sequence-connector {
margin-bottom: 5px;
}
#mismatches .form-control {
text-align: center;
height: 30px;
padding: 0;
}
"

# SHINY APP ####

runApp(shinyApp(
  ui = fluidPage(
    tags$head(tags$style(css)),
    textInput("sequence", "Enter sequence", "acggatcaag"),
    uiOutput("mismatches"),
    textOutput("oligo2")
  ),
  server = function(input, output, session) {
    output$mismatches <- renderUI({
      sequence <- toupper(input$sequence)
      #print(sequence)
      sequence= isolate_LNA(sequence)
      #print(sequence)
      
      lapply(1:length(sequence), function(x) {
        div(class = "sequence-location",
            div(sequence[x], class = "sequence-base"),
            div("|", class = "sequence-connector"),
            selectInput(paste0("base_", x), label=NULL,multiple=FALSE,selected=toupper(compgen(sequence[x])), selectize=TRUE,
                        choices= getComplementList(compgen(sequence[x]))# options=list(maxItems = 1)
            )
        )
        
      })
    })# end of mismatch input building
    

    
    output$oligo2= renderText({
      
      #check input for unexpected expressions
      oligo1=seqcheck(input$sequence)
      oligo2=NULL
      for(i in 1:nchar(oligo1)) {
        oligo2=c(oligo2,input[[paste0("base_",i)]])
      }     
      
      
      #squish, turn to lower case, and revese oligo 2 to read 5' to 3'
      oligo2=paste(oligo2,collapse="")
      oligo2=tolower(oligo2)
      oligo2=reverse_chars(oligo2)
      
    })
    
    
  }
))