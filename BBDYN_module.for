cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjmotticNAME SUBROUTINE bbdyn
cjmotticVER  versio 1.0
cjmotticPVM  last update 12.6.2018  by Juha Honkaniemi
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                    
c Simulation of bark beetle dynamics (modified from SAMBIA model by Fahse & Heurich 2011)
c BBDYN is part of the simulation framework and described in detail in:
c
c Honkaniemi, J., Ojansuu, R., Kasanen, R., & Heliövaara, K. 2018. 
c Interaction of disturbance agents on Norway spruce: a mechanistic model of bark beetle dynamics integrated in simulation framework WINDROT.
c in press Ecological Modelling.
c
c This is the main subroutine that is called by the WINDROT framework connecting the forest dynamics and disturbance models together
c This .for file contains the main subroutine bbdyn and all the other subroutines and functions needed to run it
c
c subroutine bbredistri : redistributes the bark beetles from their initial landing coordinates to nearby trees
c subroutine bb_reproduction : Colonization of trees and bark beetle reproduction in succesfully colonized trees 
c subroutine indexx1 : Indexes an array in a specified order    
c     
c function defth : Defence threshold based on tree vigour index 
c function Abark : Calculates suitable bark area for bark beetles for each tree
c function diarel : Calculates diameter of tree at height l based on stem taper curve by Laasasenaho (1982)
c function antag_curve : Probability curve for the antagonist effect
c function dheart : Calculates the heartwood diameter at diameter d
c

      subroutine bbdyn(tree,n,newtree,npuu,
     *    iseed,barkbeetle,bbphase,popbb,newpopbb)
      implicit none
	
!Subroutine inputs
	real*4 tree(50,npuu)			!Tree list including all trees  of the simulated area + buffer zone
						!           1 TreeID
						!           2 Tree sp. (1=Scots pine Pinus sylvestris, 2=Norway spruce Picea abies, 3=Silver birch Betula pendula, 4=Downy birch Betula pubescens, 5=Others)
						!           3 Tree diameter at breast height 1.3m, dbh
						!           4 x-coordinate 
						!           5 y-coordinate
						!           6 Root system radius, rr (until root diameter of 2mm)
						!           7 Living status,  -1 = living or >=0 = time from death
						!           8 Heterobasidion infection status 0 = healthy, 1 = Heterobasidion annosum, 2 = Heterobasidion parviporum
						!           9 Infected root system radius
						!           10 Time from infection, -1 = no infection or >=0 = time from infection
						!			11 Heterobasidion genet id
						!			12 Time from the infection reached stem base
                          !           13 Location (0=inside the study area, 1=buffer zone)
                          !           14 Tree type in Motti (1=live, 2=stump, 3=dead)
						!			15 Height of the decay column in stem, hdecay
						!			16 Width of the decay column in stem, ddecay
						!			17 Relative growth reduction due Heterobasidion (0-1, 1=no growth reduction) 
						!			18 Thinning method: %=0, random=1
						!			19 Tree diameter increment, cm/year
						!			20 Time from tree death
						!			21 TFlist (running numbering)
                          !           22 Tree age
                          !           23 Living crown height
                          !           24 Tree height
						!			25 -
						!			26 Decay radius in heartwood
      					!			27 Decay radius in sapwood
                          !           28 Heartwood diameter when (25) was 0
                          !           29 Radius of decay
                          !           30 Time from harvest
						!!!!!!!!!!!!!!!HWIND wind disturbance module variables
                          !           31 Root radius for HWIND 
                          !           32 Approximated decay diameter at breast height for HWIND wind model
                          !           33 Distance of the tree from stand edge of the wind direction
                          !           34 Critical wind speed for uprooting, wind speed at 10 m height in open for tree to fall down (m/s)
                          !           35 Critical wind speed for stem breakeage, wind speed at 10 m height in open for tree to break (m/s)
                          !           36 Upwind gap size
                          !           37 Wind damage status, 0=standing, 1=fallen, 2= broken
                          !           38 Maximum wind velocity for the year t 
						!!!!!!!!!!!!!!!BBRISK bark beetle module variables
                          !           39 Temporary storage of bark beetle population 
                          !           40 Number of bark beetles in the tree
                          !           41 Status of tree regarding bark beetles: 0=Healthy, 1=Failed attack, 2=Reproduction tree, 3=Dead by bark beetles
                          !           42 Colonizable bark area, Abark
                          !           43 Number of bark beetle offspring at the time step
                          !           49 Antagonist status (1=tree occupied by antagonists,0=no antagonists)
						
						
      real*4 barkbeetle(20,3000000)	!list for bark beetle individuals prior simulation step
      real*4 newtree(50,npuu)		    !Tree list including all trees  of the simulated area + buffer zone
	 
	integer npuu 					!number of trees in the tree list
	integer iseed					!seed number for random number generator
      integer n						!number of trees in the tree list				
	integer popbb					!Number of bark beetle generations from previous time step
	integer newpopbb      			!Number of bark beetle generations after the current time step
      integer bbphase                 !switch for bbdyn to determine which dispersal stage to run, 2=initialization after wind event and 1=normal stage after initialization. Determined in and provided by WINDROT framework.

!Subroutine outputs
		!Updates the tree list newtree for MOTTI and other disturbance models
      real*4 newbb(20,3000000) 		!list for bark beetle individuals after simulation step

!Helping variables
!Dispersal
      real*4 sexratio					!bark beetle population sex ratio
      integer ncohort                 !Number of beetles in a flying cohort, initial stage
      integer popul                   !Number of males in the initial beetle population popbb, bbphase=2
      integer populcoh                !Number of males per cohort x in the initial beetle population popbb, bbphase=2
	integer pop                     !Number of males in the subsequent beetle population newpopbb, bbphase=1
	integer popcoh                  !Number of males in the subsequent beetle population newpopbb, bbphase=1
      real*4 side						!study area side length (m)
	real*4 buffer  					!buffer zone width (m)
      real*4 xland                    !random x landing coordinate for beetle
	real*4 yland                    !random y landing coordinate for beetle
      real*4 xland0                   !zeroed x landing coordinate for beetle (base is the "mother tree")
	real*4 yland0                   !zeroed y landing coordinate for beetle (base is the "mother tree")
	real*4 D                        !Diffusion coefficient for bark beetle flight
	real*4 D1                       !Flying group "Strong" diffusion coefficient for bark beetle flight
	real*4 D2                       !Flying group "Average" diffusion coefficient for bark beetle flight
	real*4 D3                       !Flying group "Weak" diffusion coefficient for bark beetle flight
      real*4 tD                       !Time of dispersal
	real*4 ploc                     !Probability for beetle from cohort x to land on (x,y)
	real*4 apuloc                   !random number from equal distribution to compare against ploc
!Colonization and reproduction
      real*4 Ndef                     !Tree resistance threshold 
      real*4 Bmin                     !Minimum bark thickness for bark beetles                
	real*4 Bmax                     !Maximum bark thickness for bark beetles
	real*4 Abark                    !Colonizable bark area
      real*4 ragg                     !Aggregation catchment area radius
      real*4 rantagg                  !Anti-aggregation catchment area radius
      real*4 wintermort               !Winter mortality rate of beetles

	integer i,j,k,jj 			    !looping variables
      parameter nbb=300000            !Maximum amount of bark beetles per simulated hectare   
      
c Update the tree list to newtree to work with it in the subroutine      

      do i=1,n
		do j=1,50
			newtree(j,i)=tree(j,i)
		end do
      end do

cccccccccccccccccccccccc
ccc Model parameters ccc
cccccccccccccccccccccccc

c Population sexratio
      sexratio=0.5    	!sexratio in the bark beetle population, input is the total population, but only the dispersal of male beetles is simulated

c Optimal bark thickness minimum and maximum
      Bmin=2.5			!minimum bark thickness for bark beetle reproduction
      Bmax=999      	    !maximum bark thickness for bark beetle reproduction, 999 = no limit

c Kirjanpainajakoiraiden lentomatkat ryhmittäin  
      D1=644           	!average dispersal range of flight cohort 1 & 2 - "strong" beetles, m2 season-1
      D2=286           	!average dispersal range of fitness group 3 & 4 - "average" beetles, m2 season-1
      D3=72            	!average dispersal range of fitness group 5 & 6 - "weak" beetles, m2 season-1
      ncohort=6			!number of flight cohorts
      ragg=13.5		    !radius for aggregation catchment area, m
      rantagg=4.5		    !radius for anti-aggregation repelling area, m
      
c Threshold values for different actions related to tree's individual carrying capacity (näille arvot kirjallisuudesta)      
      wintermort=0.4    	!wintertime mortality Pollak 1975, Austarå & Midtgaard 1986     
      
c Area dimensions
      side=100 		    ! Length of effective simulation area size, m
	buffer=8  		    ! buffer
      
ccccccccccccccccccccccccc
ccc Bark beetle model ccc
ccccccccccccccccccccccccc
  
c Calculates the colonizable bark area suitable for bark beetles for each tree    
      do i=1,n
       if((newtree(14,i).eq.1.or.newtree(14,i).eq.3).and.         !trees that are alive or dead
     *         newtree(7,i).le.1.and.newtree(41,i).le.1)then      !and time from death is less than one year and are not previously succesfully colonized by bark beetles
        if(newtree(2,i).eq.2)then                                 !and are Norway spruce
        newtree(42,i)=Abark(newtree(3,i),                         !can be used attacked and used for reproduction by bark beetles
     *   newtree(24,i),Bmin,Bmax)                                 !here, the colonizable bark area is calculated for such trees
        end if
       end if
      end do                 

c Bark beetle stock population distribution in windthrown trees
      if(bbphase.eq.2) then                   !This is the first phase of the bark beetle attack, this part is run on the first year of bark beetle dynamics after the wind event
	
		popul=int(popbb*sexratio +0.5)      !we simulate only males, but assume that the catch from trap groups represents the whole population with 50:50 sexratio
		populcoh=popul*1/ncohort            !bark beetle population of a single cohort, we assume equal distribution of population into the cohorts

      do j=1,ncohort                          !looping over cohorts from 1-6
          do k=1,populcoh                     !looping over the bark beetle population in each cohort

c Calculate the landing coordinates, wind disturbance is always assumed for the lowest side of the simulation area		
          xland=100-buffer+(side+2*buffer)*ran(iseed)         !thus x-coordinates for landing vary
          yland=side                                          !but the y-coordinate stays the same, as the beetles are first spread on the warm forest edge where most of the wind-damaged trees are
      
          call bbredistri(newtree,npuu,n,ragg,                !from these landing coordinates, beetles are redistributed to trees on the forest edge
     *                    rantagg,j,xland,yland)
      
          end do
          
          !Updating the treelist on the bark beetles      
		do i=1,n
          newtree(40,i)=newtree(40,i)+newtree(43,i)            !number of beetles in the tree
          newtree(39,i)=newtree(40,i)                          !number of beetles in the tree for temporary storage                  
          newtree(43,i)=0                                      !zero the number of beetles here before using it for the next cohort and reproduction subroutine
		end do
      
      end do

c Bark beetle colonization and reproduction in the succesfully colonized trees
      call bb_reproduction(npuu,n,newtree,nbb,
     *                        barkbeetle,popbb,sexratio,iseed)
     
c Owerwinter mortality - We assume that bark beetles either overwinter in the tree or nearby the tree. 
c Here the population is summed from the parent generation and offspring and then the overwinter mortality is applied  
       do i=1,n
         if(newtree(40,i).gt.0.or.newtree(43,i).gt.0)then
         newtree(40,i)=int((newtree(40,i)+newtree(43,i))*(1-wintermort)) !over winter surviving population (previous generation and offspring) forms the next parent generation for the next time step
         end if
       end do
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Bark beetle dispersal from forest edge trees after the initial attack
      else if(bbphase.eq.1)then    !simulation is directed here by WINDROT framework wher ethe switch is changed with the wind event, first year of wind exposure is bbphase=2 and onwards bbphase=1


      do i=1,n
          newtree(43,i)=0                     !zero the storage column for bark beetle populations prior to any calculations
          newtree(39,i)=newtree(40,i)         !save the number of beetles in the temporary storage for the dispersal flight
      end do

!Simulate the flight of the bark beetles. We assume that the initial coordinates are of the tree where the beetle attacked in previous year.
      do j=1,ncohort                          !looping over cohorts from 1-6
        do i=1,n                              !looping over the tree list, where we get the i) number of beetles in the tree and ii) coordinates of the tree
            if(newtree(2,i).eq.2)then         !if tree is Norway spruce
            if(newtree(39,i).gt.0)then        !if bark beetle population is greater than 0
              pop=0                           !zero the population variable for the tree  
              popcoh=0                        !zero the cohort population variable for the tree 
              
              pop=newtree(39,i)               !population from the tree is derived from the temporary storage column
              popcoh=int(pop/ncohort+.5) 	    !population of each cohort flying off from the tree
              if(j.eq.6) popcoh=newtree(40,i) !last cohort will get the remaining amount of beetles to avoid any beetles left behind
              
            do jj=1,popcoh                    !looping over the bark beetle population in each cohort
               
              !First, determine the flight radius for the current flying cohort 
              if(j.eq.1.or.j.eq.2) D=D1  	!First 2 cohorts are "strong"
              if(j.eq.3.or.j.eq.4) D=D2   !Next 2 cohorts are "average"
              if(j.eq.5.or.j.eq.6) D=D3   !Last 2 cohorts are "weak"
              
!Randomly pick the landing coordinate location for the beetle within the simulated area + buffer zone         
145           xland=100-buffer+(side+2*buffer)*ran(iseed)
              yland=100-buffer+(side+2*buffer)*ran(iseed)
              
              !zero the coordinates based on the starting point for the flight
              xland0=abs(xland-newtree(4,i))
              yland0=abs(yland-newtree(5,i))
              
              if(xland0.eq.0.and.yland0.eq.0) go to 145       !beetle can't land on the "mother tree"
              
!Calculate the probability for the landing coordinates p(x,y) based on gaussian kernel random walk              
              tD=1                                            !time period is 1 season, D diffusion coefficient is m2 season-1
              ploc=exp(-(xland0**2+yland0**2)/(4*D*tD))       !calculation of the probability, depends on the cohort spesific diffusion coefficient D
              apuloc=0
              apuloc=ran(iseed)                               !random number from equal distribution for comparison against the probability for landing     
              if(ploc.lt.apuloc) then
                  go to 145
              end if
              
              
!If still out of the simulated area, then randomly give coordinates. This assumes that the nearby forests are similar and have equal amount of outflow of bark beetles from the area.        
              if(xland.lt.side.or.xland.gt.side*2) then
                  xland=100-buffer+(side+2*buffer)*ran(iseed)
                  yland=100-buffer+(side+2*buffer)*ran(iseed)
              end if
              
              if(yland.lt.side.or.yland.gt.side*2) then
                  xland=100-buffer+(side+2*buffer)*ran(iseed)
                  yland=100-buffer+(side+2*buffer)*ran(iseed)
              end if
			  
!Redistribution of bark beetles from their original landing coordinates to trees due to bark beetle pheromone aggregation and anti-aggregation   
			call bbredistri(newtree,npuu,n,ragg,
     *                            rantagg,j,xland,yland)
              
			end do

! After the cohort is redistributed, then reduce the number of beetles in the cohort from the bark beetle population of the "mother tree"
			newtree(40,i)=newtree(40,i)-popcoh
              
			end if
			end if
        end do
        
!Updating the treelist on the bark beetles  
		do i=1,n
          newtree(40,i)=newtree(40,i)+newtree(43,i)   !number of beetles in the tree
          newtree(43,i)=0                             !zero the number of beetles here before using it for the next cohort and reproduction subroutine
		end do
      
      end do
             
c And finally, outside the whole loop to avoid interference, let's update the temporary storage again for bark beetle population
      do i=1,n
          newtree(39,i)=0                         !zero it
          newtree(39,i)=newtree(40,i)             !and update
      end do
                 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Bark beetle colonization of new trees and mortality due to unsuccesful colonization
       call bb_reproduction(npuu,n,newtree,nbb,
     *                        newbb,newpopbb,sexratio,iseed)
          
c Owerwinter mortality - We assume that bark beetles either overwinter in the tree or nearby the tree. 
c Here the population is summed from the parent generation and offspring and then the overwinter mortality is applied  
       do i=1,n
         if(newtree(40,i).gt.0.or.newtree(43,i).gt.0)then
         newtree(40,i)=int((newtree(40,i)+newtree(43,i))*(1-wintermort))   !over winter surviving population (previous generation and offspring) forms the next parent generation for the next time step
         end if
       end do
   
      end if !END OF THE BBPHASE IF CLAUSE
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccc          
      return   !Return to MOTTI 
      end
         
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjmottic SUBROUTINE bbredistri
cjmottic	Redistributes the bark beetles from their initial landing coordinates to nearby trees
c			due to aggregation and anti-aggregation pheromones emitted by attacking bark beetles
c			See Honkaniemi et al 2018 "Submodels (i)"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

      subroutine bbredistri(newtree,npuu,n,ragg,
     * rantagg,j,xland,yland)
      
	implicit none 
      
!Subroutine inputs
	real*4 	newtree(50,npuu),       !Tree list including all trees  of the simulated area + buffer zone
     * 			ragg,			    !Model parameter for pheromone aggregation area radius ragg
     * 			rantagg,		    !Model parameter for pheromone anti-aggregation area radius rantagg
     *			xland,				!Initial x coordinate for bark beetle individual landing
     *			yland				!Initial y coordinate for bark beetle individual landing
	 
	integer npuu,n, j				!npuu=number of trees, n=number of trees, j=flying cohorts of bark beetles (1-2=strong, 3-4=average and 5-6=weak)
      
!Subroutine outputs
	!Updated tree list newtree regarding the number of beetles attacking the tree
      
!Helping variables
	real*4  agglist(8,npuu)       	!aggregation list, list of trees  
      real*4  distbb					!distance from bark beetle initial landing coordinates to trees within the aggreation distance
      integer ind2(npuu)			    !indexed list of trees from nearest to furthest around the beetle landing coordinates
      real*4  aggrlim, antaggrlim	    !number of beetles in a tree needed to i) start aggregating other beetles and ii) to repell other beetles with anti-aggregation pheromones
      real*4  eta					    !distance between anti-aggregation tree and other trees - beetle can't move to trees within the anti-aggregation area 			
      integer defth					!Defence threshold Ndef = Number of bark beetles needed to overcome tree resistance
	integer ii,k,kk,l				!loop indices
      
!Re-distribution of bark beetles based on aggregation and anti-aggregation phereomones or in their absence to the nearest tree              

! Define first the trees within the maximum final movement distance from the landing spot
          kk=0  
          do ii=1,n                                               !looping through the whole tree list of the stand
              if(newtree(2,ii).eq.2)then                          !if tree is Norway spruce
              if((newtree(14,ii).eq.1.or.newtree(14,ii).eq.3).and.!if tree is alive or dead for maximum one year and that it hasn't been succesfully colonized by bark beetles before
     *         newtree(7,ii).le.1.and.newtree(41,ii).le.1)then   
           distbb=sqrt((newtree(4,ii)-xland)**2+     				!distance from the random landing coordinates of the bark beetle individual to the tree
     *        (newtree(5,ii)-yland)**2)
           
           if(distbb.lt.ragg)then 							!if distance to the tree is less than the aggregation distance, then the tree is listed to lie within the aggregation area for bark beetle
           kk=kk+1
           agglist(1,kk)=ii					                    !tree number
           agglist(2,kk)=distbb				                    !distance to the tree
           agglist(3,kk)=newtree(40,ii)		                    !Number of bark beetles in the tree				
           agglist(4,kk)=newtree(4,ii)		                    !x-coordinate of the tree
           agglist(5,kk)=newtree(5,ii)		                    !y-coordinate of the tree
           agglist(6,kk)=newtree(42,ii)		                    !colonizable bark area of the tree
           agglist(7,kk)=newtree(17,ii) 	                        !possible growth reduction in the tree due to root rot
           agglist(8,kk)=newtree(37,ii) 	                        !Wind damage status, 0=standing, 1=fallen, 2= broken
           end if
              end if
              end if
          end do
              
       if(kk.gt.0)then  !if such trees are found then the bark beetle is arranged to these as follows:
      
!First the trees within the aggregation distance are indexed from closest to furthest from the original bark beetle landing coordinates
          do ii=1,kk
           call indexx1(kk,agglist(2,1:kk),ind2(1:kk))   
          end do
              
!Aggregation and anti-aggregation of beetles to trees
		
! Anti-aggregation of flying cohorts STRONG, strong beetles are not affected by the aggregation pheromones       
		if(j.eq.1.or.j.eq.2)then                                        !only flying cohorts 1 and 2
            do k=1,kk                   
				antaggrlim=defth(newtree(3,agglist(1,ind2(k))),newtree(19,agglist(1,
     *     		ind2(k))),newtree(2,agglist(1,ind2(k))),
     *     		newtree(22,agglist(1,ind2(k))))                         !defence threshold of the tre. When that is reached, bark beetles start to produce anti-aggregating pheromones 
 
				if(agglist(3,ind2(k)).ge.antaggrlim)then		        !if the nearest tree is within anti-aggregation range
				   do l=1,kk
					eta=sqrt((agglist(4,ind2(k))-agglist(4,ind2(l)))**2+
     *        		(agglist(5,ind2(k))-agglist(5,ind2(l)))**2)         !distance between the anti-aggregating tree and other trees in the vicinity. Beetle can't be moved to a tree within the anti-aggregation range
					if(eta.gt.rantagg.and.
     *            	agglist(3,ind2(l)).lt.antaggrlim)then               !tree either not anti-aggregation tree or outside the range of other anti-aggregating trees
					newtree(43,agglist(1,ind2(l)))=newtree(43,agglist(1,ind2(l)))+1 !then add a beetle to the number of bark beetles in the tree
					go to 212     							            !return the updated tree list
					end if
				   end do
				end if               
            end do   
		end if
      
! Aggregation of AVERAGE and WEAK flying cohorts, anti-aggregation does not affect them
		if(j.gt.2)then !flying cohorts 3-6
            do k=1,kk                  
				aggrlim=1 				                                !the first beetle attacking the tree starts to produce aggregation pheromone. As the STRONG beetles are not affected, the number after those cohorts is often higher than that
				if(agglist(3,ind2(k)).ge.aggrlim)then                   !if number of beetles in the tree bigger than the limit
				newtree(43,agglist(1,ind2(k)))=newtree(43,agglist(1,ind2(k)))+1 !then add a beetle to the number of bark beetles in the tree
                go to 212                                                 !return the updated tree list
				end if
            end do
		end if     
              
!Last, if there are no aggregating or anti-aggregating trees, the beetle is moved to the nearest tree          
			do k=1,kk        	
          aggrlim=1 			                                            !the first beetle attacking the tree starts to produce aggregation pheromone. As the STRONG beetles are not affected, the number after those cohorts is often higher than that
				if(agglist(3,ind2(k)).lt.aggrlim.and.
     *        		agglist(3,ind2(k)).ge.0)then   						!if beetles less than the aggregation limit has attacked the tree
					newtree(43,agglist(1,ind2(k)))=newtree(43,agglist(1,ind2(k)))+1 	!then add a beetle to the number of bark beetles in the tree
					go to 212 		                                    !return the updated tree list
				end if               
            end do  
      end if
         
212   return  !return the updated tree list
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjmottic SUBROUTINE bb_reproduction
cjmottic	Colonization of trees and bark beetle reproduction in succesfully colonized trees
c			See Honkaniemi et al 2018 "Submodels (ii & iii)"
c
c			For the antagonists see also:
c			Fahse and Heurich 2011 
c			Honkaniemi et al 2018 "Submodels (iv)"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bb_reproduction(npuu,n,newtree,nbb,
     * barkbeetle,popbb,sexratio,iseed)
      
      implicit none
!Subroutine input
      integer npuu,n,nbb,popbb  
      real*4
     * newtree(50,npuu),			        !Tree list including all trees  of the simulated area + buffer zone
     * barkbeetle(20,3000000),            !list for bark beetle individuals prior simulation step
     * sexratio						    !Bark beetle population sexratio
!Subroutine outputs
	integer newpopbb				    !New number of bark beetles within the subroutine
      !Also updated tree list newtree regarding the number of beetles in the tree
	 
!Helping variables
	integer Nfem					    !Number of females per tree
      integer Nf1					        !Number of offspring per tree
      integer popbb_temp			        !Temporary number of bark beetles within the subroutine	
	integer i,j					        !looping
      integer defenceth				    !Defence threshold Ndef = Number of bark beetles needed to overcome tree resistance
	integer defth					    !Function for defence threshold
      integer iseed   				    !seed number for random number generator
	real*4  Dfem					    !Density of females per 100cm2 of bark
	
!Antagonists
	real*4 Nanttrees				    !Total number of trees within the antagonist area from tree i
	real*4 Nbbkills				        !Number of bark beetle killed trees in the previous year within the antagonist area
	real*4 pfind                        !Probability for the antagonists to find bark beetle reproduction trees
	real*4 Pbbkills					    !Proportion of infested trees in the antagonist area
      real*4 apuant					    !random number from equal distribution (0,1) to compare against  pfind
	real*4 rpfind				        !Antagonist area where they can still find the bark beetle reproduction trees
	real*4 antag_eff				    !Efficiency of the antagonists to kill larvae
	real*4 antag_curve			        !Function for probability curve for the antagonist effect
	real*4 antag_dist				    !Distance of trees between the antagonist area
	real*4 pi						    !pi
	real*4 p0						    !lower bound of antagonist curve
	real*4 p1						    !upper bound of antagonist curve

!Setting up the variables      
      pi=3.1415926
      p0=0.03            			        !lower bound of antagonist curve
      p1=0.999            			    !upper bound of antagonist curve
      antag_eff=0.95       			    !efficiency of the antagonists to kill larvae

   
! First the colonization of the trees if the defence threshold is exceeded during the dispersal    
      do i=1,n
        if(newtree(2,i).eq.2)then       		    !If tree is Norway spruce
         if(newtree(40,i).gt.0)then				!If the number of attacking bark beetles is more than 0
          defenceth=defth(newtree(3,i),newtree(19,i),newtree(2,i),
     *     newtree(22,i))             			!Calculation of the defence threshold for tree     
          if(newtree(40,i).gt.defenceth) then 	!If defence threshold overcome
          newtree(41,i)=2  						!Change of trees infestation status
          end if
         end if
        end if
      end do
      
! Next step is the reproduction only in succesfully colonized trees, HOWEVER

! Before that we want to kill the beetles in trees which were not succesfully colonized
      
      do i=1,n
        if(newtree(2,i).eq.2)then                     !If tree is Norway spruce
         if(newtree(40,i).gt.0.and.newtree(41,i).le.1)then      !if tree has beetles, but has not been succesfully colonized 
          newtree(40,i)=0  						    !Kill the beetles from those trees
          newtree(41,i)=1  						    !Change the tree status to attacked, but failed. Considered as healthy, but helps accounting the dynamics
         end if
        end if
      end do
          
! Then is turn for the beetle reproduction 
      do i=1,n
        if(newtree(2,i).eq.2)then
         if(newtree(41,i).eq.2.)then
          Nfem=newtree(40,i)*2  			            !number of beetles per tree is the number of males, thus to get the number of females, we assume that every male has two females Nmale*2=Nfem
          Dfem=Nfem/newtree(42,i)*100		            !Density of females per 100cm2 to match the function from Andebrandt et al 1985
          Nf1=Dfem*38.8*exp(-0.87*(Dfem**0.45))*(newtree(42,i)/100)   !Number of offspring produced per tree, see Andebrandt et al 1985
          newtree(43,i)=int(Nf1*sexratio)             !Number of offspring times sexratio (0.5 by default) to keep track on only males
         end if
        end if
      end do
 
! N.B. Intraspecific competition is taken into account in the reproduction function, winter mortality taken into account later in the end of the simulation step in subroutine bbdyn

! Here, last step is to take into account the effect of antagonists on larval survival
!     See Fahse & Heurich (2011) for further details
   
! First, calculate the area where antagonists can in theory find bark beetle reproduction trees
      rpfind=sqrt(2500/pi)                     !That area is 0.25 ha = 2500 m2

! Second, calculate the number of trees within this antagonist area as well as trees that have been killed by bark beetles in the previous year. This increases then the probability for antagonists to find new reproduction trees.
      do i=1,n							        !looping over the whole tree list
        if(newtree(2,i).eq.2)then		            !If tree i is Norway spruce
          if(newtree(43,i).gt.0)then   	        !If tree i has been used for reproduction Nf1>0
            Nbbkills=0 					        !zero the number of trees killed by bark beetles in the previous years within the antagonist dispersal area
            Nanttrees=0					        !zero the total number of trees within the antagonist dispersal area
          do j=1,n						        !inner loop over the whole tree list
            if(newtree(2,j).eq.2)then	            !If tree j is Norway spruce
            if(newtree(7,j).le.1)then  	        !If tree j is alive or dead for less than a year
            antag_dist=sqrt((newtree(4,i)-newtree(4,j))**2+
     *        (newtree(5,i)-newtree(5,j))**2)		!Calculate the distance between tree i and j
            if(antag_dist.lt.rpfind)then 		!If tree j is within the antagonist area from the tree i (current reproduction tree) 
            Nanttrees=Nanttrees+1					!Add one to the number of trees within the antagonist area			
            if(newtree(41,j).ge.2.and.newtree(7,j).le.1)then     !Trees within antagonist area that were killed by bark beetles in the previous year
            Nbbkills=Nbbkills+1
            end if
            end if
            end if
            end if
          end do

! Calculate the proportion of infested trees in the vicinity, this will be an input variable to calculate the find probability
           if(Nanttrees.gt.0)then		        !If trees in the antagonist area
            Pbbkills=Nbbkills/Nanttrees	
           else
            Pbbkills=0					    !if no trees, probability is 0
           end if
          
		  pfind=antag_curve(p1,p0,Pbbkills)	!Calculate the find probability according to Fahse & Heurich 2011	
          
            apuant=ran(iseed)				    !random number from equal distribution (0,1)
           if(apuant.lt.pfind) then			!compare the random number to the calculated finding probability
            newtree(49,i)=1					!if passed, the tree i is occupied by antagonists
            newtree(43,i)=int(newtree(43,i)*(1-antag_eff)) !then also the offspring is killed by the antagonists, with a survival probability of (1-antag_eff) (NOTE! in this study the survival was 0.05)
           else
            newtree(49,i)=0					!if tree is not occupied by antagonists, then nothing happens either to the bark beetle offspring
           end if
          end if
        end if
      end do  
   
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjmottic FUNCTION defth
cjmottic	Defence threshold based on tree vigour index 
c			See Honkaniemi et al (2018) "Submodels (ii) Colonization and mortality of trees"
c			See Mulock and Christiansen (1986) for further details for original study
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function defth(dbh,id13,pl,age)
      implicit none
!Function input
	real*4 	dbh,		!Tree diameter at breast height
     *		id13,		!Annual tree growth at breast height
     *		pl,			!tree species code (1=Scots pine, 2=Norway spruce)	
     *		age			!Tree age
!Function output
	integer defth 		!Defence threshold Ndef = Number of bark beetles needed to overcome tree resistance
!Helping variables
	real*4 	pi, 		!pi
     *		BAi,		!Basal area of annual growth
     *		SA,			!Sapwood area
     *		dhea,       !Heart wood diameter
     *		dheart      !function for calculating the heartwood diameter

      pi=3.1415926
      dhea=dheart(pl,dbh,age)         !Calculating the heartwood diameter for the tree
               
      BAi=pi*dbh**2-pi*(dbh-id13)**2  !Calculate the basal area of annual growth
      SA=pi*dbh**2-pi*dhea**2         !Calculate the sapwood area 
      
      defth=90*exp(0.24*(BAi/SA*100)) !Tree vigour index tvig=BAi/SA transformed to percents as in Mulock & Christiansen 1986    
      defth=defth*(dbh/20)            !Mulock and christiansen use dbh=20cm as a base value for their equation, so transform here to real dbh value from the Mulock & Christiansen 1986 equation
      
      return
      end
      
	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjmottic FUNCTION Abark
cjmottic	Calculates suitable bark area for bark beetles for each tree
c			See Honkaniemi et al. 2018 "Submodels (iii)"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      function Abark(dbh,h,Bmin,Bmax)
	implicit none
!Function input
      real*4 dbh, 			!tree diameter at breast height
     *		h, 			    !tree height
     *     	Bmin, 			!minimum bark thickness for bark beetle reproduction
     *		Bmax			!maximum bark thickness for bark beetle reproduction
	  
!Function output
	real*4 Abark			!colonizable bark area for bark beetle reproduction
	  
!Helping variables
	real*4 relx,  		    !Relative tree height
     * 	 lBmin, 		    !Height at Bmin 
     *	 lBmax, 		    !Height at Bmax
     * 	 Blj,			    !Bark thickness at height within tree l
     * 	 btI			    !rounded bark thickness for comparison with minimum and maximum values
	  
	integer l  			    !Height within tree, here used in a loop to calculate the relative tree height
      integer i				!looping variable 
      
	real*4 pi				!pi
	real*4 relxmin,relxmax  !relative height values within loop, temporary storage 
	  
      real*4 l1,l2  		    !heights of upper and lower bounds of bark area
	real*4 rl1,rl2		    !diameters at upper and lower bounds of bark area
	real*4 p1,p2 			!perimeters at upper and lower bounds of bark area
      
	real*4 hAbark			!Height of the colonizable bark area
	real*4 h10			    !Height of the 1/10 of the colonizable bark area
	real*4 len			    !Height of the 1/10 of the colonizable bark area
	real*4 A10			    !1/10 of the colonizable bark area
	real*4 A				!Cumulative colonizable bark area within loop		
      real*4 lBmax1           !initial value within loop for the height of maximum bark thickness
      
	real*4 diarel   		!function for diameter of tree at height l based on stem taper curve by Laasasenaho (1982)
      
!Setting parameters within function
      pi=3.1415926
      A=0.
      lBmin=0.
      lBmax=999.
      lBmax1=0.
      Abark=0.
      
!Divide the stem to 1000 parts and calculate the bark thicknesses in each 0.1% relative height of the stem
      do l=0,1000
      relx=(1-l/1000.)
      
      Blj=(19.97778*relx+0.19083*(dbh*relx)-40.69507*relx**2+13.17111*
     * relx**3+41.36845*relx**5-54.48194*relx**8+26.91165*relx**13)    !see Honkaniemi et al 2018 Eq 5
        
      btI=anint(Blj*10)/10  	        !Rounding the bark thickness for comparison against min and max parameter values
      
      if(btI.eq.Bmin) then  	        !If bark thickness at l equals minimum parameter value Bmin then the relative point of height is recorded
          lBmin=h*(1-relx)
          relxmin=relx
      end if
      if(btI.eq.Bmax) then 		    !If bark thickness at l equals maximum parameter value Bmax then the relative point of height is recorded
          lBmax1=h*(1-relx)
          if(lBmax1.lt.lBmax)then 
           lBmax=h*(1-relx)
           relxmax=relx
          end if
      end if
      end do
      
      if(lBmax.le.0..or.lBmax.eq.999) lBmax=0.    !making the value 0 if too low or too high values for some reason
      if(lBmin.le.0..or.lBmin.eq.999) lBmin=0.	!making the value 0 if too low or too high values for some reason
      
795      hAbark=(h-(h-lBmin)-lBmax) 		        !calculates the height of the bark area between 
      
!Calculates the colonizable bark area in 10 equal length trapezoid parts: See Honkaniemi et al. 2018 Eq 7
      do i=1,10				
      len=0.
      len=0.1*hAbark
      l1=0.
      l2=0.
      
      l1=(lBmax+(i-1)*len)
      l2=(lBmax+i*len)
      h10=l2-l1   
        
      rl1=diarel(dbh,h,l1,2.)
      rl2=diarel(dbh,h,l2,2.)

      p1=pi*rl1 					    !perimeter at lower bound of the 10th bark area
      p2=pi*rl2 					    !perimeter at lower bound of the 10th bark area
      
      A10=(p1+p2)/2*h10*100			!10th of the colonizable bark area calculated as an area of trapezoid
      A=A+A10
      A10=0.
      end do
      
      Abark=A   					    !function returns the value of colonizable bark area for bark beetle reproduction in single tree
      
	return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjmottic FUNCTION diarel
cjmottic	Calculates diameter of tree at height l based on stem taper curve by Laasasenaho
c			See Honkaniemi et al. 2018 "Submodels (iii)"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function diarel(dbh,h,l,pl)
      implicit none
!Function inputs
	real*4  dbh,			!tree diameter at breast height
     *	 	h,				!tree height
     * 	    l,				!tree height at height l
     *		pl				!tree species code (1=Scots pine, 2=Norway spruce)	
	  
!Function output
      real*4 diarel			!diameter of tree at height l 
	  
!Laasesenaho (1982) model parameters for Norway spruce
	real*4 b1,b2,b3,b4,b5,b6,b7,b8
     
!Helping variables
	real*4 x13,x  		    !relative height at dbh and relative height at l
       
      x13=1-1.3/h
      x=1-l/h

      if(pl.eq.2)then    	    !model only for Norway spruce included here as Ips typographus attacks only Norway spruce
          b1=2.3366
          b2=-3.2684
          b3=3.6513
          b4=-2.2608
          b5=0.0
          b6=2.1501
          b7=-2.7412
          b8=1.8876
      end if
      
!Model returns the relative height at dbh and relative height at l
      diarel=(dbh/(b1*x13+b2*x13**2+b3*x13**3+b4*x13**5+b5*x13**8+
     * b6*x13**13+b7*x13**21+b8*x13**34))*(b1*x+b2*x**2+b3*x**3+b4*x**5+
     * b5*x**8+b6*x**13+b7*x**21+b8*x**34) 

      return 
      end
  
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjmottic FUNCTION antag_curve
cjmottic	Probability curve for the antagonist effect
c			See Fahse and Heurich 2011 
c			See Honkaniemi et al 2018 "Submodels (iv)"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function antag_curve(pfind1,pfind0,a)
	implicit none
!Function inputs
	real*4 	pfind1,			!Probability for the trees to be found by antagonists when all trees in the vicinity are infested or dead due to bark beetles
     *		pfind0,			!Probability for the trees to be found by antagonists when none of the trees in the vicinity are infested or dead due to bark beetles
     *		a				!The proportion of infested trees in the vicinity
!Function output
	real*4  antag_curve     !Antagonism probability for the bark beetle reproduction process

!Helping variables
	real*4  c,rm,fi         !

      c=pfind0/(1-pfind0)
      rm=log(pfind1/(c*(1-pfind1)))
      fi=c*exp(a*rm)
      antag_curve=fi/(1+fi)   !Return the antagonistic probability
 
	return
      end    

cccccccccccccccccccccccccccccccccccccccccccccccc
cjmotticNIMI FUNCTION dheart
cjmottic	Calculates the heartwood diameter at diameter d
c         See Honkaniemi et al. 2014. Hmodel, a Heterobasidion annosum model for even-aged Norway spruce stands.
c             Canadian Journal of Forest Research 44:796-809. dx.doi.org/10.1139/cjfr-2014-0011
cccccccccccccccccccccccccccccccccccccccccccccccc
c
	function dheart(pl,d,age)
	implicit none
	real*4 dheart,			    ! Hearwood diameter at diameter d
     *	pl,					    ! Tree sp. (1=Scots pine Pinus sylvestris, 2=Norway spruce Picea abies, 3=Silver birch Betula pendula, 4=Downy birch Betula pubescens, 5=Others)
     *	d,					    ! Diameter of the tree above bark
     *	age					    ! Tree age
	real*4 dmm                  ! diameter in mm with bark
      real*4 dbmm                 ! diameter in mm without bark
c
	dheart=0.
      dmm=d*10
      if(pl.eq.1) then                                !equation for Scots pine
          dbmm=dmm-exp(3.4445+0.00515*dmm-0.03598*65)		! tree diameter without bark
          dheart=(-15.4+0.1580*dbmm*log(age))/10          ! heartwood diameter
      else if (pl.eq.2) then                          !equation for Scots pine
          dbmm=dmm-exp(0.6644+0.00091*dmm+0.2907*log(dmm))! tree diameter without bark
          dheart=(-15.6+0.2149*dbmm*log(age)              ! heartwood diameter
     *			-0.00124*dbmm*(log(age))**3)/10
      end if

	if(dheart.lt.0.) dheart=0.01                !if heartwood in some case smaller than 0, then leave always some heartwood in all trees
	return
      end
	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjmottic SUBROUTINE indexx1
cjmottic Indexes an array in a specified order
cjmottic Original source code from book Numerical recipes in Fortran 77 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
      subroutine indexx1(n,arr,indx)
      implicit none
      integer n,indx(n),M,NSTACK
      real*4 arr(n)
      parameter (M=7,NSTACK=50)
      !Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that arr(indx(j))
      !is in ascending order for j = 1;2; : : :;N. The input quantities n and arr are not changed.
      integer i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      real*4 a
      do j=1,n
          indx(j)=j
      enddo 
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
          do j=l+1,ir
              indxt=indx(j)
              a=arr(indxt)
              do i=j-1,l,-1
                  if(arr(indx(i)).le.a)goto 2
                  indx(i+1)=indx(i)
              enddo
              i=l-1
2             indx(i+1)=indxt
          enddo
          if(jstack.eq.0)return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
      else
          k=(l+ir)/2
          itemp=indx(k)
          indx(k)=indx(l+1)
          indx(l+1)=itemp
          if(arr(indx(l)).gt.arr(indx(ir)))then
              itemp=indx(l)
              indx(l)=indx(ir)
              indx(ir)=itemp
          endif
          if(arr(indx(l+1)).gt.arr(indx(ir)))then
              itemp=indx(l+1)
              indx(l+1)=indx(ir)
              indx(ir)=itemp
          endif
          if(arr(indx(l)).gt.arr(indx(l+1)))then
              itemp=indx(l)
              indx(l)=indx(l+1)
              indx(l+1)=itemp
          endif
          i=l+1
          j=ir
          indxt=indx(l+1)
          a=arr(indxt)
3         continue
          i=i+1
          if(arr(indx(i)).lt.a)goto 3
4         continue
          j=j-1
          if(arr(indx(j)).gt.a)goto 4
          if(j.lt.i)goto 5
          itemp=indx(i)
          indx(i)=indx(j)
          indx(j)=itemp
          goto 3
5         indx(l+1)=indx(j)
          indx(j)=indxt
          jstack=jstack+2
          if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
          if(ir-i+1.ge.j-l)then
              istack(jstack)=ir
              istack(jstack-1)=i
              ir=j-1
          else
              istack(jstack)=j-1
              istack(jstack-1)=l
              l=i
          endif
      endif
      goto 1
      end subroutine 