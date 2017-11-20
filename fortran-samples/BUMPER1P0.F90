MODULE diatom
      IMPLICIT NONE
      TYPE species_data
      REAL, DIMENSION(0:1000)::abundance=0.0 !Abundance in each training set site. abundance(0) is dummy value used in RMSE.
      REAL, DIMENSION(1000)::core_abundance !Abundance in each core sample
      REAL:: wpeak_abundance,woptimum,wtolerance,wpower,wprob_pa !probability weighted species response curve
      REAL:: WAopt,PAtol !optimum and tolerance (WA optimum, PA tolerance) 
      REAL:: mean_TS_abundance !Mean abundance given presence in training set, used for adundance prior
      REAL::prob_pres !scaler for probability of presence (multiply by prob_pa for probability of presence at optimum)
      REAL, DIMENSION(0:3,0:9,0:3,0:3,0:3):: probsrc=0.0!Probability of all SRCs.
      REAL::prob_score=0.0 !measure of fit of most likely SRC - used to evaluate alternative distributions
      REAL::max_prob_src !probability of most probable SRC
      REAL::p_optimum
      INTEGER::no_counts !Number of lakes which contain the species
      INTEGER::no_srcs !Number of significant SRCs (probability >10% most probable SRC, max_prob_src)
      CHARACTER*20:: species !species name
      END TYPE species_data
END MODULE diatom

MODULE training_set
      USE diatom
      IMPLICIT NONE
      TYPE training_data_set
      TYPE (species_data), DIMENSION(1000)::species !Data for each species in training set
      REAL, DIMENSION(1000) :: env_var !Environmental variable in each training set lake
      REAL, DIMENSION(1000,100) :: rec_env_var_pdf=0 !posterior for reconstructed env_variable in training set
      REAL, DIMENSION(1000,100) :: core_rec_env_var_pdf=0 !posterior for reconstructed env_variable in training set
      REAL, DIMENSION(1000) :: rec_env_var=0!Reconstructed env_variable in training set
      REAL, DIMENSION(1000) :: core_rec_env_var=0!Reconstructed env_variable in training set
      REAL::SRC_opt_max,SRC_opt_min,SRC_tol_max,SRC_tol_min,SRC_power_min,SRC_power_max,SRC_prob_pa_max,SRC_prob_pa_min! LIMITS OF SRCs
      REAL::SRC_abund_mult_max,SRC_abund_mult_min!! What abundance range tested? Multipliers to mean abundance in TS
      REAL::max_env_var,min_env_var!LIMITS OF TS
      REAL::rec_range_max,rec_range_min!LIMITS OF RECONSTRUCTION
      REAL::pa_component!percentage of likelihood function defined by presence absence term
      REAL::required_abundance!minimum abundance for inclusion in reconstruction
      INTEGER, DIMENSION(1000) :: richness !number of species with (abundance > required abunance ) in training set
      INTEGER, DIMENSION(1000) :: core_richness !number of species with (abundance > required abunance ) in training set
      INTEGER::SRC_i,SRC_j,SRC_k,SRC_l,SRC_m!DIMENSIONS OF SRCs TESTED (ABUNDANCE, OPTIMUM, TOLERANCE, POWER, PROB)  
      INTEGER:: no_species !total number of species
      INTEGER:: no_lakes !total number of lakes
      INTEGER:: no_core_samples !total number of core samples
      CHARACTER*20:: name
      CHARACTER*20, DIMENSION(1000):: lake_name
      CHARACTER*20, DIMENSION(1000):: core_name
      END TYPE training_data_set
END MODULE training_set

PROGRAM Bayes_Reconstruction
      USE training_set
      IMPLICIT NONE
      TYPE (training_data_set)::data_set
      INTEGER::n,i,j,k,l,m
      INTEGER::switch_core,switch_jack,switch_lake,switch_app

      PRINT*,'DATA SET NAME?'
      READ(*,*) data_set%name
      PRINT*,'CALCULATE APPARENT RMSE? 1=YES'
      READ(*,*) switch_app
      PRINT*,'CALCULATE JACKNIFED RMSEP? 1=YES'
      READ(*,*) switch_jack
      IF(switch_jack.ne.1) THEN
       PRINT*,'CALCULATE INDIVIDUAL JACKNIFED LAKE? 0=NONE >0 LAKE ID'
       READ(*,*) switch_lake
      ELSE
       switch_lake=0
      ENDIF
      PRINT*,'RECONSTRUCT CORE? 1=YES'
      READ(*,*) switch_core

      CALL SETUP(data_set)

      IF(switch_app.eq.1.or.switch_core.eq.1.or.(switch_app+switch_jack+switch_lake+switch_core).eq.0) THEN
        PRINT*,'CALCULATING RESPONSE CURVES'
        CALL CALCULATE_RESPONSE_CURVES(0,data_set)! CALCULATE SRCS FOR FULL DATA SET
        CALL OUTPUT_SRC(data_set)
      ENDIF


      IF(switch_jack.eq.1.or.switch_lake.gt.0) then
        OPEN(UNIT=3,FILE='ts_likelihood_jack.data',ACCESS='APPEND')
        WRITE (UNIT=3,FMT=200) 0,0,0,0,(data_set%rec_range_min+real(n)*(data_set%rec_range_max-&
          data_set%rec_range_min)/99.0,n=0,99)
  200   FORMAT(I4,I25,I4,I25,1P100E11.3)
        CLOSE(UNIT=3)
       ENDIF

      IF(switch_app.EQ.1) THEN
        PRINT*,'CALCULATING APPARENT RMSE'
        OPEN(UNIT=3,FILE='ts_likelihood_app.data',ACCESS='APPEND')
        WRITE (UNIT=3,FMT=200) 0,0,0,0,(data_set%rec_range_min+real(n)*(data_set%rec_range_max-&
          data_set%rec_range_min)/99.0,n=0,99)
        CLOSE(UNIT=3)
        DO i=1,data_set%no_lakes
          CALL RECONSTRUCT(1,i,data_set) !RECONSTRUCT LAKES FOR APPARENT RMSE
        ENDDO
        CALL DIAGNOSE_RMSE(data_set,0)
        CALL OUTPUT_TS_POST(data_set,0)
      ENDIF

      IF(switch_core.EQ.1) THEN
        PRINT*,'RECONSTRUCTING CORE'
        OPEN(UNIT=3,FILE='core_likelihood.data',ACCESS='APPEND')
        WRITE (UNIT=3,FMT=200) 0,0,0,0,(data_set%rec_range_min+real(n)*(data_set%rec_range_max-&
          data_set%rec_range_min)/99.0,n=0,99)
        CLOSE(UNIT=3)
        DO i=1,data_set%no_core_samples
          CALL RECONSTRUCT(0,i,data_set) !RECONSTRUCT CORE
        ENDDO
        CALL DIAGNOSE_CORE(data_set)
        CALL OUTPUT_CORE_POST(data_set)
      ENDIF

      IF(switch_jack.EQ.1) THEN
        PRINT*,'CALCULATING RMSEP'
        DO i=1,data_set%no_lakes
          CALL CALCULATE_RESPONSE_CURVES(i,data_set)! CALCULATE JACKNIFED SRCS    
          CALL RECONSTRUCT(2,i,data_set) !RECONSTRUCT JACKKNIFED LAKE
        ENDDO
        CALL DIAGNOSE_RMSE(data_set,1)
        CALL OUTPUT_TS_POST(data_set,1)
      ENDIF

      IF(switch_lake.GT.0) THEN
        CALL CALCULATE_RESPONSE_CURVES(switch_lake,data_set)! CALCULATE INIDIVIDUAL JACKNIFED SRC    
        CALL OUTPUT_SRC(data_set)
        CALL RECONSTRUCT(2,switch_lake,data_set) !RECONSTRUCT INDIVIDUAL JACKKNIFED LAKE
      ENDIF

END PROGRAM Bayes_Reconstruction


SUBROUTINE SETUP(data_set)
      USE training_set
      IMPLICIT NONE
      TYPE (training_data_set)::data_set
      INTEGER:: no_env_var,which_env_var
      INTEGER::i,j,n
      INTEGER:: ithresh,imodel
      REAL, DIMENSION(20):: read_variables
      REAL atol,ntol !for estimation of average tolerence, used in priors
      REAL sampling_density
      character*20 data_string

!TRAINING SET ENVIRONMENT
      OPEN (UNIT=3,FILE=trim(data_set%name)//'.lakes')
      READ(3,*) data_set%no_lakes
      DO i=1,data_set%no_lakes
        READ(3,*) data_set%lake_name(i),data_set%env_var(i)
      ENDDO
      CLOSE(UNIT=3)

!SPECIES LIST 
      OPEN (UNIT=3,FILE=trim(data_set%name)//'.species')
      READ(3,*)data_set%no_species
      DO j=1,data_set%no_species
        READ(3,*) data_set%species(j)%species
      ENDDO
      CLOSE(UNIT=3)

!TRAINING SET %COUNTS        
      OPEN (UNIT=3,FILE=trim(data_set%name)//'.ts.counts')
      DO i=1,data_set%no_lakes
        READ(3,*) (data_set%species(j)%abundance(i),j=1,data_set%no_species)
      ENDDO
      CLOSE(UNIT=3)

!INPUT CORE DATA
      OPEN (UNIT=3,FILE=trim(data_set%name)//'.core.counts')
      READ(3,*)data_set%no_core_samples
      DO i=1,data_set%no_core_samples
        READ(3,*) data_set%core_name(i),&
                  (data_set%species(j)%core_abundance(i),j=1,data_set%no_species)
      ENDDO
      CLOSE (UNIT=3)

!DUMMY LAKE 0 - SRC(0) CALL INCLUDES ALL SPECIES
      DO j=1,data_set%no_species
        data_set%species(j)%abundance(0)=100. !just a logical requirement so all species are considered in SRC calculations
      ENDDO

!CALCULATE MEAN ABUNDANCE (GIVEN PRESENCE) OF EACH SPECIES WITHIN TRAINING SET AND NUMBER OF OCCURRANCES.
!THESE USE THE ENTIRE TRAINING SET, INCLUDING LAKES WITH USE_LAKE==0
      sampling_density=0.0
      DO i=1,data_set%no_species
        data_set%species(i)%mean_TS_abundance=sum(data_set%species(i)%abundance(1:data_set%no_lakes))
        data_set%species(i)%no_counts=0
        DO j=1,data_set%no_lakes
          if(data_set%species(i)%abundance(j)>0) THEN
            data_set%species(i)%no_counts=1+data_set%species(i)%no_counts
            sampling_density=sampling_density+1.0
          ENDIF            
        ENDDO
        IF(data_set%species(i)%no_counts.eq.0) THEN
          data_set%species(i)%mean_TS_abundance=0.0
        ELSE     
          data_set%species(i)%mean_TS_abundance= &
            data_set%species(i)%mean_TS_abundance/data_set%species(i)%no_counts
        ENDIF
      ENDDO
!      print*,
!      print*,"Sampling density",sampling_density/data_set%no_species
!      print*,"Mean richness",sampling_density/data_set%no_lakes

!CALCULATE WA OPTIMA AND TOLERANCES TO INFORM PRIORS
      DO i=1,data_set%no_species
        data_set%species(i)%WAopt=0.0
        IF(data_set%species(i)%no_counts.gt.0) THEN !need 1 count for an estimate of optimum
          DO j=1,data_set%no_lakes
            data_set%species(i)%WAopt=data_set%species(i)%WAopt+ &
             data_set%species(i)%abundance(j)*data_set%env_var(j)
          ENDDO
          data_set%species(i)%WAopt=data_set%species(i)%WAopt/ &
           sum(data_set%species(i)%abundance(1:data_set%no_lakes))
        ENDIF
      ENDDO
! TOLERANCE BASED ONLY ON PRESENCE
      atol=0.0
      ntol=0.0
      DO i=1,data_set%no_species
        data_set%species(i)%PAtol=0.0
        if(data_set%species(i)%no_counts.gt.1) then !need 2 counts for an estimate of tolerance. Increase this?
          ntol=ntol+1.
          DO j=1,data_set%no_lakes
            if(data_set%species(i)%abundance(j).gt.0.0) then
              data_set%species(i)%PAtol=data_set%species(i)%PAtol+ &
               (data_set%env_var(j)-data_set%species(i)%WAopt)**2
            endif
          ENDDO
          data_set%species(i)%PAtol=sqrt(data_set%species(i)%PAtol/ &
           data_set%species(i)%no_counts)
          atol=atol+data_set%species(i)%PAtol
        endif
      ENDDO
      atol=atol/ntol
!      print*,"indicative tolerance",atol

!CALCULATE ENVIRONMENTAL RANGE
! THIS USES ALL LAKES, INCLUDING THOSE WITH USE_LAKE==0
      data_set%max_env_var=maxval(data_set%env_var(1:data_set%no_lakes))
      data_set%min_env_var=minval(data_set%env_var(1:data_set%no_lakes))

!ARRAY SIZE OF PROBSRC - CHANGE ARRAY SIZE OF PROBSRC IF THESE CHANGE
      data_set%SRC_i=3  !abundance
      data_set%SRC_j=9 !optimum
      data_set%SRC_k=3  !tolerance
      data_set%SRC_l=3  !power
      data_set%SRC_m=3  !prob_pa

      print*,
      print*,"PRESENCE/ABSENCE MODEL (YES=1)"
      read*,imodel
      data_set%pa_component=0.5
      if(imodel.eq.1) data_set%pa_component=1.0
      print*,
      print*,"APPLY 2% THRESHOLD? (YES=1)"
      read*,ithresh
      data_set%required_abundance=0.001
      if(ithresh.eq.1) data_set%required_abundance=2.0

!box prior for environmental variable
      data_set%rec_range_min=data_set%min_env_var-atol*6.0  
      data_set%rec_range_max=data_set%max_env_var+atol*6.0

!box priors for SRC
      data_set%SRC_opt_max=data_set%max_env_var+atol
      data_set%SRC_opt_min=data_set%min_env_var-atol
      data_set%SRC_tol_max=atol*3.0
      data_set%SRC_tol_min=atol*2.0/3.0
      data_set%SRC_power_max=1.0
      data_set%SRC_power_min=0.2
      data_set%SRC_prob_pa_max=2.5
      data_set%SRC_prob_pa_min=1.0  
      data_set%SRC_abund_mult_max=2.5
      data_set%SRC_abund_mult_min=1.0

!PROBABILITY OF PRESENCE, USED TO DEFINE prob_pa prior - maximum=0.99/prob_pa_max to cap prior at 0.99
      DO i=1,data_set%no_species
        data_set%species(i)%prob_pres=&
          min((real(data_set%species(i)%no_counts)/real(data_set%no_lakes)),&
          (0.99/data_set%SRC_prob_pa_max))
      ENDDO

!sample richness
      data_set%richness(:)=0
      do i=1,data_set%no_lakes
        do j=1,data_set%no_species
          if(data_set%species(j)%abundance(i).gt.data_set%required_abundance) then
            data_set%richness(i)=data_set%richness(i)+1
          endif
        enddo
      enddo
!sample richness core
      data_set%core_richness(:)=0
      do i=1,data_set%no_core_samples
        do j=1,data_set%no_species
          if(data_set%species(j)%core_abundance(i).gt.data_set%required_abundance) then
            data_set%core_richness(i)=data_set%core_richness(i)+1
          endif
        enddo
      enddo

!INITIALISE OUTPUT FILES
      OPEN(UNIT=3,FILE='output.data',STATUS='REPLACE')
      CLOSE(UNIT=3)
      OPEN(UNIT=3,FILE='weighted_src.data',STATUS='REPLACE')
      WRITE(3,*)' ID   OPTIMUM TOLERANCE ABUNDANCE     POWER PROB_PRES   NO_SRCS  NO_LAKES'
      CLOSE(UNIT=3)
      OPEN(UNIT=3,FILE='reconstruct_ts_app.data',STATUS='REPLACE')
      CLOSE(UNIT=3)
      OPEN(UNIT=3,FILE='reconstruct_ts_jack.data',STATUS='REPLACE')
      CLOSE(UNIT=3)
      OPEN(UNIT=3,FILE='reconstruct_core.data',STATUS='REPLACE')
      CLOSE(UNIT=3)
      OPEN(UNIT=3,FILE='ts_likelihood_app.data',STATUS='REPLACE')
      CLOSE(UNIT=3)
      OPEN(UNIT=3,FILE='ts_likelihood_jack.data',STATUS='REPLACE')
      CLOSE(UNIT=3)
      OPEN(UNIT=3,FILE='core_likelihood.data',STATUS='REPLACE')
      CLOSE(UNIT=3)
    
END SUBROUTINE SETUP


SUBROUTINE CALCULATE_RESPONSE_CURVES(jackknife,data_set)
!if jackknife=0, all TS lakes included - for apparent RMSE and core reconstructions
      USE training_set
      IMPLICIT NONE
      TYPE (training_data_set)::data_set    
      REAL:: prob
      REAL:: Log_probability
      REAL:: max_abundance,min_abundance,test_abundance
      REAL:: max_optimum,min_optimum,test_optimum
      REAL:: max_tolerance,min_tolerance,test_tolerance
      REAL:: max_power,min_power,test_power
      REAL:: max_prob_pa,min_prob_pa,test_prob_pa 
      REAL::max_prob_src
      REAL:: TS_abundance,abundance
      REAL::model_fit
      REAL::Gaussian
      INTEGER:: i,j,k,l,m,n,o
      INTEGER::jackknife

      max_optimum=data_set%SRC_opt_max
      min_optimum=data_set%SRC_opt_min
      max_tolerance=data_set%SRC_tol_max
      min_tolerance=data_set%SRC_tol_min
      max_power=data_set%SRC_power_max
      min_power=data_set%SRC_power_min
      max_prob_pa=data_set%SRC_prob_pa_max
      min_prob_pa=data_set%SRC_prob_pa_min

      DO n=1,data_set%no_species
        data_set%species(n)%probsrc=epsilon(prob)
        data_set%species(n)%wpeak_abundance=0
        data_set%species(n)%woptimum=0
        data_set%species(n)%wtolerance=0
        data_set%species(n)%wpower=0
        data_set%species(n)%wprob_pa=0
      ENDDO

      DO n=1,data_set%no_species
!dont bother calculating if less than required abundance in jacknifed lake.
!note jacknife=0 for non-jacknife calculations (dummy abundance=100)
        if(data_set%species(n)%abundance(jackknife)<data_set%required_abundance) cycle
        max_abundance=data_set%species(n)%mean_TS_abundance*data_set%SRC_abund_mult_max
        min_abundance=data_set%species(n)%mean_TS_abundance*data_set%SRC_abund_mult_min
!LOOP SRCS
        DO i=0,data_set%SRC_i !LOOP TO TEST PEAK ABUNDANCE
          test_abundance=min_abundance+REAL(i)*(max_abundance-min_abundance)/REAL(data_set%SRC_i)
          DO j=0,data_set%SRC_j ! LOOP TO TEST OPTIMIMUM
            test_optimum=min_optimum+REAL(j)*(max_optimum-min_optimum)/REAL(data_set%SRC_j)
            DO k=0,data_set%SRC_k !LOOP TO TEST TOLERANCE
              test_tolerance=min_tolerance+REAL(k)*(max_tolerance-min_tolerance)/REAL(data_set%SRC_k)
              DO l=0,data_set%SRC_l !LOOP TO TEST TOLERANCE
                test_power=min_power+REAL(l)*(max_power-min_power)/REAL(data_set%SRC_l)                    
                DO m=0,data_set%SRC_m !LOOP TO TEST TOLERANCE
                  test_prob_pa=min_prob_pa+REAL(m)*(max_prob_pa-min_prob_pa)/REAL(data_set%SRC_m)

                  do o=1,data_set%no_lakes
                    if(o==jackknife) cycle !JACK-KNIFED LAKE NOT IN SRC CALCULATION
                    TS_abundance=data_set%species(n)%abundance(o)
                    data_set%species(n)%probsrc(i,j,k,l,m)=&
                      data_set%species(n)%probsrc(i,j,k,l,m)+&
                      Log_probability(data_set%env_var(o),test_abundance,test_optimum,&
                      test_tolerance,test_power,test_prob_pa,TS_abundance,&
                      data_set%species(n)%prob_pres)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        data_set%species(n)%probsrc=data_set%species(n)%probsrc-maxval(data_set%species(n)%probsrc)
        data_set%species(n)%probsrc=EXP(data_set%species(n)%probsrc)
        data_set%species(n)%probsrc=data_set%species(n)%probsrc/sum(data_set%species(n)%probsrc)
      ENDDO



      DO n=1,data_set%no_species
        max_prob_src=0.0
        data_set%species(n)%no_srcs=0
        DO i=0,data_set%SRC_i
          DO j=0,data_set%SRC_j
            DO k=0,data_set%SRC_k
              DO l=0,data_set%SRC_l
                DO m=0,data_set%SRC_m                        
! MOST PROBABLE SRC
                  IF(data_set%species(n)%probsrc(i,j,k,l,m)>max_prob_src) THEN
                    max_prob_src=data_set%species(n)%probsrc(i,j,k,l,m)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        data_set%species(n)%max_prob_src=max_prob_src
!NUMBER OF SIGNIFICANT SRCSs
        DO i=0,data_set%SRC_i
          DO j=0,data_set%SRC_j
            DO k=0,data_set%SRC_k
              DO l=0,data_set%SRC_l
                DO m=0,data_set%SRC_m
                  if(data_set%species(n)%probsrc(i,j,k,l,m)>max_prob_src*0.01) &
                    data_set%species(n)%no_srcs=data_set%species(n)%no_srcs+1
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO        
      ENDDO

!Calculate prob weighted SRCs (output data, not used in calculation)
      DO n=1,data_set%no_species
        if(data_set%species(n)%abundance(jackknife)<data_set%required_abundance) cycle !no need to calculate SRC if species<required abundance
        max_abundance=data_set%species(n)%mean_TS_abundance*data_set%SRC_abund_mult_max
        min_abundance=data_set%species(n)%mean_TS_abundance*data_set%SRC_abund_mult_min
        DO i=0,data_set%SRC_i !LOOP TO TEST PEAK ABUNDANCE
          test_abundance=min_abundance+REAL(i)*(max_abundance-min_abundance)/REAL(data_set%SRC_i)            
          DO j=0,data_set%SRC_j ! LOOP TO TEST OPTIMIMUM
            test_optimum=min_optimum+REAL(j)*(max_optimum-min_optimum)/REAL(data_set%SRC_j)
            DO k=0,data_set%SRC_k !LOOP TO TEST TOLERANCE
              test_tolerance=min_tolerance+REAL(k)*(max_tolerance-min_tolerance)/REAL(data_set%SRC_k)
              DO l=0,data_set%SRC_l !LOOP TO TEST TOLERANCE
                test_power=min_power+REAL(l)*(max_power-min_power)/REAL(data_set%SRC_l)                    
                DO m=0,data_set%SRC_m !LOOP TO TEST TOLERANCE
                  test_prob_pa=min_prob_pa+REAL(m)*(max_prob_pa-min_prob_pa)/REAL(data_set%SRC_m)                    
                  data_set%species(n)%wpeak_abundance=data_set%species(n)%wpeak_abundance+&
                    test_abundance*data_set%species(n)%probsrc(i,j,k,l,m)
                  data_set%species(n)%woptimum=data_set%species(n)%woptimum+&
                    test_optimum*data_set%species(n)%probsrc(i,j,k,l,m)
                  data_set%species(n)%wtolerance=data_set%species(n)%wtolerance+&
                    test_tolerance*data_set%species(n)%probsrc(i,j,k,l,m)
                  data_set%species(n)%wpower=data_set%species(n)%wpower+&
                    test_power*data_set%species(n)%probsrc(i,j,k,l,m)
                  data_set%species(n)%wprob_pa=data_set%species(n)%wprob_pa+&
                    test_prob_pa*data_set%species(n)%probsrc(i,j,k,l,m)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

END SUBROUTINE CALCULATE_RESPONSE_CURVES


SUBROUTINE RECONSTRUCT(ID,which,data_set)
! ID=0 reconstructs core sample
! ID=1 reconstructs TS lake (all lakes) 
! ID=2 reconstructs TS lake (jacknifed)
! NB NO DIFFERENCE BETWEEN 1&2 EXCEPT THEY OUTPUT TO DIFFERENT FILES
      USE training_set
      IMPLICIT NONE
      TYPE (training_data_set)::data_set    
      REAL:: Log_probability
      REAL:: max_env_var
      REAL:: min_env_var
      REAL:: test_env_var
      REAL::max_abundance,min_abundance
      REAL::prob(100)
      REAL,DIMENSION(1000)::maximum,pdf_total !for normalisation
      REAL,DIMENSION(1000)::maximum_likelihood,minimum_likelihood,total_likelihood !for likelihood normalisation
      REAL:: peak,opt,tol,power,prob_pa,log_p,tol_p
      REAL::likelihood_fn(1000,100) !LIKELIHOOD FUNCTION FOR LAKE
      REAL::likelihood_fn_y(1000,100),likelihood_fn_p(1000,100)!Two likelihood function terms
      REAL::env_var_interval
      REAL::abundance,sample_abundance !tested expected fractional abundance and actual fractional abundance
      REAL::Gaussian
      INTEGER::ID
      INTEGER::which 
      INTEGER:: i,j,k,l,m,n,o,loop

      max_env_var=data_set%rec_range_max
      min_env_var=data_set%rec_range_min
      env_var_interval=(max_env_var-min_env_var)/99.
                
      likelihood_fn=epsilon(likelihood_fn)
      likelihood_fn_y=epsilon(likelihood_fn_y)
      likelihood_fn_p=epsilon(likelihood_fn_p)

! LIKELIHOOD COUNT TERM
      DO n=1,data_set%no_species

        IF(ID.GT.0) THEN 
          sample_abundance=data_set%species(n)%abundance(which)
        ELSE
          sample_abundance=data_set%species(n)%core_abundance(which)
        ENDIF            

        IF(sample_abundance<data_set%required_abundance) CYCLE

        max_abundance=data_set%species(n)%mean_TS_abundance*data_set%SRC_abund_mult_max
        min_abundance=data_set%species(n)%mean_TS_abundance*data_set%SRC_abund_mult_min

          DO i=0,data_set%SRC_i
            peak=min_abundance+REAL(i)*(max_abundance-min_abundance)/REAL(data_set%SRC_i)            
            DO j=0,data_set%SRC_j
              opt=data_set%SRC_opt_min+REAL(j)*(data_set%SRC_opt_max-data_set%SRC_opt_min)/&
                REAL(data_set%SRC_j)
              DO k=0,data_set%SRC_k
                tol=data_set%SRC_tol_min+REAL(k)*(data_set%SRC_tol_max-&
                  data_set%SRC_tol_min)/REAL(data_set%SRC_k)
                DO l=0,data_set%SRC_l
                  power=data_set%SRC_power_min+REAL(l)*(data_set%SRC_power_max-&
                    data_set%SRC_power_min)/REAL(data_set%SRC_l)
                  DO m=0,data_set%SRC_m
!only SRCs with >1% of max prob are considered for speed
                    if(data_set%species(n)%probsrc(i,j,k,l,m)<data_set%species(n)%max_prob_src*0.01) CYCLE
                    prob_pa=data_set%SRC_prob_pa_min+&
                      REAL(m)*(data_set%SRC_prob_pa_max-data_set%SRC_prob_pa_min)/&
                      REAL(data_set%SRC_m)
                    test_env_var=min_env_var
                    DO o=1,100
                      log_p=Log_probability(test_env_var,peak,opt,tol,power,&
                        prob_pa,sample_abundance,data_set%species(n)%prob_pres)
                      likelihood_fn_y(n,o)=likelihood_fn_y(n,o)+exp(log_p)*data_set%species(n)%probsrc(i,j,k,l,m)
                       test_env_var=test_env_var+env_var_interval
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
      ENDDO

      maximum_likelihood=MAXVAL(likelihood_fn_y,DIM=2)
      minimum_likelihood=MINVAL(likelihood_fn_y,DIM=2)    
      total_likelihood=SUM(likelihood_fn_y,DIM=2)

      DO n=1,data_set%no_species
        DO o=1,100
          likelihood_fn_y(n,o)=likelihood_fn_y(n,o)/total_likelihood(n)
        ENDDO
      ENDDO

!LIKELIHOOD PRESENCE TERM

      DO n=1,data_set%no_species
        IF(ID.GT.0) THEN 
          sample_abundance=data_set%species(n)%abundance(which)
        ELSE
          sample_abundance=data_set%species(n)%core_abundance(which)
        ENDIF            
        if(sample_abundance<data_set%required_abundance) cycle
!!!???
        if(sample_abundance<0.0001) cycle! no likelihood term if species not present (as likelihood is pa term anyway)
        test_env_var=min_env_var
        DO o=1,100
          DO i=0,data_set%SRC_i
            DO j=0,data_set%SRC_j
              opt=data_set%SRC_opt_min+REAL(j)*(data_set%SRC_opt_max-data_set%SRC_opt_min)/REAL(data_set%SRC_j)
              DO k=0,data_set%SRC_k
                tol=data_set%SRC_tol_min+REAL(k)*(data_set%SRC_tol_max-data_set%SRC_tol_min)/REAL(data_set%SRC_k)
                DO l=0,data_set%SRC_l  
                  power=data_set%SRC_power_min+REAL(l)*(data_set%SRC_power_max-&
                    data_set%SRC_power_min)/REAL(data_set%SRC_l)
                  tol_p=tol/sqrt(power)
                  DO m=0,data_set%SRC_m
                    if(data_set%species(n)%probsrc(i,j,k,l,m)<data_set%species(n)%max_prob_src*0.01) CYCLE
                    prob_pa=data_set%SRC_prob_pa_min+REAL(m)*(data_set%SRC_prob_pa_max-data_set%SRC_prob_pa_min)/&
                      REAL(data_set%SRC_m)
                    prob_pa=prob_pa*data_set%species(n)%prob_pres
                    likelihood_fn_p(n,o)=likelihood_fn_p(n,o)+(prob_pa*&
                      exp(-((test_env_var-opt)*(test_env_var-opt)/(2*tol_p*tol_p)))*&
                      data_set%species(n)%probsrc(i,j,k,l,m))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          test_env_var=test_env_var+env_var_interval
        ENDDO
      ENDDO
      total_likelihood=SUM(likelihood_fn_p,DIM=2)

!NORMALISE LIKELIHOOD PRESENCE FUNCTION

      DO n=1,data_set%no_species
        DO o=1,100
          likelihood_fn_p(n,o)=likelihood_fn_p(n,o)/total_likelihood(n)
        ENDDO
      ENDDO

!COMBINE PRESENCE AND COUNT LIKELIHOOD FUNCTIONS

      DO n=1,data_set%no_species
        DO o=1,100
          likelihood_fn(n,o)=(likelihood_fn_p(n,o)*data_set%pa_component)+&
            (likelihood_fn_y(n,o)*(1.0-data_set%pa_component))
        ENDDO
      ENDDO

!OUTPUT LIKELIHOOD FUNCTIONS
  300 FORMAT(I4,A25,I4,A25,100E11.4)

      IF(ID.EQ.0) THEN
        OPEN (UNIT=3,FILE='core_likelihood.data',ACCESS='APPEND')
        DO n=1,data_set%no_species
          IF(minimum_likelihood(n)<maximum_likelihood(n)) THEN
            WRITE(UNIT=3,FMT=300) which,trim(data_set%core_name(which)),n,trim(data_set%species(n)%species),&
                                  (likelihood_fn(n,o),o=1,100)
          ENDIF
        ENDDO
        CLOSE(UNIT=3)
      ENDIF

      IF(ID.EQ.1) THEN
        OPEN (UNIT=3,FILE='ts_likelihood_app.data',ACCESS='APPEND')
        DO n=1,data_set%no_species
          IF(minimum_likelihood(n)<maximum_likelihood(n)) THEN
          WRITE(UNIT=3,FMT=300) which,trim(data_set%lake_name(which)),n,trim(data_set%species(n)%species),&
                                (likelihood_fn(n,o),o=1,100)
          ENDIF
        ENDDO
        CLOSE(UNIT=3)
      ENDIF

      IF(ID.EQ.2) THEN
        OPEN (UNIT=3,FILE='ts_likelihood_jack.data',ACCESS='APPEND')
        DO n=1,data_set%no_species
          IF(minimum_likelihood(n)<maximum_likelihood(n)) THEN
          WRITE(UNIT=3,FMT=300) which,trim(data_set%lake_name(which)),n,trim(data_set%species(n)%species),&
                                (likelihood_fn(n,o),o=1,100)
          ENDIF
        ENDDO
        CLOSE(UNIT=3)
      ENDIF


!COMBINE SPECIES LIKELIHOOD FUNCTIONS AND NORMALISE
      IF(ID.GT.0) THEN
        data_set%rec_env_var(which)=0.0
        DO o=1,100
          data_set%rec_env_var_pdf(which,o)=0.0
          DO n=1,data_set%no_species
            sample_abundance=data_set%species(n)%abundance(which)
            IF(sample_abundance<data_set%required_abundance) cycle
            data_set%rec_env_var_pdf(which,o)=data_set%rec_env_var_pdf(which,o)+log(likelihood_fn(n,o))
          ENDDO
    	ENDDO
        maximum=MAXVAL(data_set%rec_env_var_pdf,DIM=2)
        DO o=1,100
          data_set%rec_env_var_pdf(which,o)=data_set%rec_env_var_pdf(which,o)-maximum(which)
          data_set%rec_env_var_pdf(which,o)=exp(data_set%rec_env_var_pdf(which,o))
        ENDDO
        pdf_total=sum(data_set%rec_env_var_pdf,DIM=2)
        test_env_var=min_env_var
        DO o=1,100
          data_set%rec_env_var(which)=data_set%rec_env_var(which)+test_env_var*&
          data_set%rec_env_var_pdf(which,o)/pdf_total(which)
          test_env_var=test_env_var+env_var_interval
        ENDDO
	print*,which,data_set%env_var(which),data_set%rec_env_var(which),trim(data_set%lake_name(which))
      ELSE
        data_set%core_rec_env_var(which)=0.0
        DO o=1,100
          data_set%core_rec_env_var_pdf(which,o)=0.0
          DO n=1,data_set%no_species
            sample_abundance=data_set%species(n)%core_abundance(which)
            IF(sample_abundance<data_set%required_abundance) cycle
    	    data_set%core_rec_env_var_pdf(which,o)=data_set%core_rec_env_var_pdf(which,o)+log(likelihood_fn(n,o))
          ENDDO
        ENDDO
        maximum=MAXVAL(data_set%core_rec_env_var_pdf,DIM=2)
        DO o=1,100
          data_set%core_rec_env_var_pdf(which,o)=data_set%core_rec_env_var_pdf(which,o)-maximum(which)
          data_set%core_rec_env_var_pdf(which,o)=exp(data_set%core_rec_env_var_pdf(which,o))
        ENDDO
        pdf_total=sum(data_set%core_rec_env_var_pdf,DIM=2)
        test_env_var=min_env_var
        DO o=1,100
          data_set%core_rec_env_var(which)=data_set%core_rec_env_var(which)+test_env_var*&
          data_set%core_rec_env_var_pdf(which,o)/pdf_total(which)
          test_env_var=test_env_var+env_var_interval
        ENDDO
        print*,which,data_set%core_rec_env_var(which)
      ENDIF
      
END SUBROUTINE RECONSTRUCT


SUBROUTINE OUTPUT_SRC(data_set)
      USE training_set
      IMPLICIT NONE    
      TYPE (training_data_set)::data_set
      INTEGER:: n
  100 FORMAT(I4,5F10.3,2I10,a20)
      OPEN (UNIT=3,FILE='weighted_src.data',ACCESS='APPEND')
      DO n=1,data_set%no_species
        WRITE (UNIT=3,FMT=100) n,data_set%species(n)%woptimum,data_set%species(n)%wtolerance,&
          data_set%species(n)%wpeak_abundance,data_set%species(n)%wpower,&
          data_set%species(n)%wprob_pa*data_set%species(n)%prob_pres,&
          data_set%species(n)%no_srcs,data_set%species(n)%no_counts,trim(data_set%species(n)%species)
      ENDDO
      CLOSE (UNIT=3)
END SUBROUTINE OUTPUT_SRC

SUBROUTINE OUTPUT_TS_POST(data_set,switch) !0=APP, 1=JACK
      USE training_set
      IMPLICIT NONE    
      TYPE (training_data_set)::data_set
      INTEGER:: n,m,switch
  200 FORMAT(I4,A25,1P100E11.3)
      if(switch==0) then
        OPEN (UNIT=3,FILE='reconstruct_ts_app.data',ACCESS='APPEND')
      else
        OPEN (UNIT=3,FILE='reconstruct_ts_jack.data',ACCESS='APPEND')
      endif
      WRITE (UNIT=3,FMT=200) 0,'SITE',(data_set%rec_range_min+real(n)*(data_set%rec_range_max-&
      data_set%rec_range_min)/99.0,n=0,99)
      DO m=1,data_set%no_lakes
        WRITE (UNIT=3,FMT=200) m,trim(data_set%lake_name(m)),(data_set%rec_env_var_pdf(m,n),n=1,100)
      ENDDO
      CLOSE (UNIT=3)
END SUBROUTINE OUTPUT_TS_POST

SUBROUTINE OUTPUT_CORE_POST(data_set)
      USE training_set
      IMPLICIT NONE    
      TYPE (training_data_set)::data_set
      INTEGER:: n,m
  200 FORMAT(I4,A25,1P100E11.3)
      OPEN (UNIT=3,FILE='reconstruct_core.data',ACCESS='APPEND')
      WRITE (UNIT=3,FMT=200) 0,'SAMPLE',(data_set%rec_range_min+real(n)*(data_set%rec_range_max-&
        data_set%rec_range_min)/99.0,n=0,99)
      DO m=1,data_set%no_core_samples
        WRITE (UNIT=3,FMT=200) m,trim(data_set%core_name(m)),&
                               (data_set%core_rec_env_var_pdf(m,n),n=1,100)
      ENDDO
      CLOSE (UNIT=3)
END SUBROUTINE OUTPUT_CORE_POST

REAL FUNCTION Gaussian(variable,max,opt,tol)
      REAL variable,max,opt,tol
      Gaussian=max*exp(-(variable-opt)*(variable-opt)/(2*tol*tol))
END FUNCTION Gaussian


REAL FUNCTION Log_probability(env_var,SRC_peak,SRC_opt,SRC_tol,SRC_power,SRC_prob_pa,counted,prob_presence)
      IMPLICIT NONE
      REAL::Gaussian
      REAL:: env_var,SRC_peak,SRC_opt,SRC_tol,SRC_power,SRC_prob_pa,counted,expected
      REAL::p_max
      REAL::prob_presence
!PROBABILITY OF PRESENCE
      p_max=SRC_prob_pa*prob_presence
      IF(counted<0.001) THEN
        Log_probability=LOG(1-(p_max*exp(-(env_var-SRC_opt)*(env_var-SRC_opt)*SRC_power/(2*SRC_tol*SRC_tol))))
        RETURN
      ELSE
        expected=Gaussian(env_var,SRC_peak,SRC_opt,SRC_tol)
        IF(expected<0.0001) THEN
          Log_probability=Log(epsilon(Log_probability))
        RETURN
        ENDIF
        Log_probability=LOG(p_max)-(env_var-SRC_opt)*(env_var-SRC_opt)*SRC_power/(2*SRC_tol*SRC_tol)
      ENDIF
! PROBABILITY OF COUNT 
! EXPONENTIAL
      Log_probability=Log_probability-LOG(expected*(1-exp(-100/expected)))-counted/expected
END FUNCTION Log_probability

SUBROUTINE DIAGNOSE_RMSE(data_set,switch)
!switch=0 app, switch=1 jack, identical calulations, flag only to idenitify output
      USE training_set
      IMPLICIT NONE
      TYPE (training_data_set)::data_set
      REAL,DIMENSION(100)::posterior
      REAL:: max_env_var
      REAL:: min_env_var
      REAL::env_var_interval
      REAL:: test_env_var
      REAL::sum_x_2,sum_x,sum_y,sum_xy
      REAL::LLSESP,average_bias,intercept,RMSE
      REAL::post_width,sum_pdf,mean_post_width
      REAL::R2,meanlake,meanrec,Exy,sx,sy
      REAL::mbias(10),mbiascount(10),thresh
      INTEGER::lake,o,biascount
      INTEGER switch

      OPEN(UNIT=99,FILE='output.data',ACCESS='APPEND')

      mbias=0.0
      mbiascount=0.0

      max_env_var=data_set%rec_range_max
      min_env_var=data_set%rec_range_min
      env_var_interval=(max_env_var-min_env_var)/99.0
!for least squares fit
      sum_x=0.0
      sum_x_2=0.0
      sum_y=0.0
      sum_xy=0.0
      DO lake=1,data_set%no_lakes
        sum_x=sum_x+data_set%env_var(lake)
        sum_x_2=sum_x_2+data_set%env_var(lake)*data_set%env_var(lake)
      ENDDO
      average_bias=0.0
      LLSESP=1.0
      RMSE=0.0
      DO lake=1,data_set%no_lakes
        data_set%rec_env_var(lake)=0.0
        test_env_var=min_env_var
        DO o=1,100
          posterior(o)=data_set%rec_env_var_pdf(lake,o) !no prior
          data_set%rec_env_var(lake)=data_set%rec_env_var(lake)+test_env_var*posterior(o)
          test_env_var=test_env_var+env_var_interval
        ENDDO
        data_set%rec_env_var(lake)=data_set%rec_env_var(lake)/sum(posterior)
!LEAST SQUARES FIT
        sum_y=sum_y+data_set%rec_env_var(lake)
        sum_xy=sum_xy+data_set%env_var(lake)*data_set%rec_env_var(lake)
        RMSE=RMSE+(data_set%rec_env_var(lake)-data_set%env_var(lake))**2
!MAX BIAS
        thresh=data_set%min_env_var+(data_set%max_env_var-data_set%min_env_var)*0.1
        do biascount=1,10
          if(data_set%env_var(lake).le.thresh) then
            mbias(biascount)=mbias(biascount)+data_set%rec_env_var(lake)-data_set%env_var(lake)
            mbiascount(biascount)=mbiascount(biascount)+1
            exit
          endif
          thresh=thresh+(data_set%max_env_var-data_set%min_env_var)*0.1
        enddo
      ENDDO  
      do biascount=1,10
        mbias(biascount)=mbias(biascount)/mbiascount(biascount)    
      enddo
!R2
      meanlake=sum_x/float(data_set%no_lakes)
      meanrec=sum_y/float(data_set%no_lakes)
      Exy=0.0
      sx=0.0
      sy=0.0
      do lake=1,data_set%no_lakes
        Exy=Exy+(data_set%env_var(lake)-meanlake)*(data_set%rec_env_var(lake)-meanrec)
        sx=sx+(data_set%env_var(lake)-meanlake)**2
        sy=sy+(data_set%rec_env_var(lake)-meanrec)**2
      enddo
      R2=(Exy/(sqrt(sx)*sqrt(sy)))**2
!solve gradient of least squares fit (measured v reconstructed)
      LLSESP=((data_set%no_lakes*sum_xy)-(sum_x*sum_y))/((data_set%no_lakes*sum_x_2)-(sum_x*sum_x))
      intercept=(sum_y-LLSESP*sum_x)/data_set%no_lakes
      average_bias=(sum_y-sum_x)/data_set%no_lakes
      RMSE=sqrt(RMSE/data_set%no_lakes)       

      IF(switch==1) then
        print*,'JACKNIFED STATISTICS'
        print*,'RMSEP',RMSE
      else
        print*,'FITTED MODEL STATISTICS'
        print*,'RMSE',RMSE
      endif
      print*,'R2',R2
      print*,'Average bias',average_bias
      print*,'Maximum bias',maxval(mbias)
      print*,'LLSESP',LLSESP-1.0
      print*,'Intercept',Intercept

      IF(switch==1) then
        write(99,*) 'JACKNIFED STATISTICS'
        write(99,*) 'RMSEP',RMSE
      else
        write(99,*) 'FITTED MODEL STATISTICS'
        write(99,*) 'RMSE',RMSE
      endif
      WRITE(99,*)'R2',R2
      WRITE(99,*)'Average bias',average_bias
      WRITE(99,*)'Maximum bias',maxval(mbias)
      WRITE(99,*)'LLSESP',LLSESP-1.0
      WRITE(99,*)'Intercept',Intercept

      mean_post_width=0.0
      WRITE(99,*) 
      WRITE(99,*) '                         LAKE  MEASURED     RECON POST_WDTH  RICHNESS'  
      DO lake=1,data_set%no_lakes
        test_env_var=min_env_var
        post_width=0.0
        sum_pdf=0.0
        DO o=1,100
          post_width=post_width+data_set%rec_env_var_pdf(lake,o)*(test_env_var-data_set%rec_env_var(lake))**2
          sum_pdf=sum_pdf+data_set%rec_env_var_pdf(lake,o)
          test_env_var=test_env_var+env_var_interval
        ENDDO
        post_width=sqrt(post_width/sum_pdf)
        mean_post_width=mean_post_width+post_width
        WRITE(99,990) trim(data_set%lake_name(lake)),data_set%env_var(lake),data_set%rec_env_var(lake),&
                      post_width,data_set%richness(lake)
      ENDDO
      mean_post_width=mean_post_width/data_set%no_lakes
      WRITE(99,*)
      WRITE(99,*) 'MEAN POSTERIOR WIDTH',mean_post_width
      print*,'MEAN POSTERIOR WIDTH',mean_post_width
  990 FORMAT(A25,3F10.3,I10)
      CLOSE(UNIT=99)
END SUBROUTINE DIAGNOSE_RMSE

SUBROUTINE DIAGNOSE_CORE(data_set)
      USE training_set
      IMPLICIT NONE
      TYPE (training_data_set)::data_set
      REAL:: max_env_var
      REAL:: min_env_var
      REAL::env_var_interval
      REAL:: test_env_var
      REAL::sum_pdf,post_width
      INTEGER::core,o

      OPEN(UNIT=99,ACCESS='APPEND',FILE='output.data')

      max_env_var=data_set%rec_range_max
      min_env_var=data_set%rec_range_min
      env_var_interval=(max_env_var-min_env_var)/99.0

      DO core=1,data_set%no_core_samples
        data_set%core_rec_env_var(core)=0.0
        test_env_var=min_env_var
        DO o=1,100
          data_set%core_rec_env_var(core)=data_set%core_rec_env_var(core)+test_env_var*data_set%core_rec_env_var_pdf(core,o) !no prior
          test_env_var=test_env_var+env_var_interval
        ENDDO
      ENDDO      

      DO core=1,data_set%no_core_samples
        sum_pdf=0.0
        DO o=1,100
          sum_pdf=sum_pdf+data_set%core_rec_env_var_pdf(core,o)
        ENDDO
        data_set%core_rec_env_var(core)=data_set%core_rec_env_var(core)/sum_pdf
      ENDDO

      WRITE(99,*) 'CORE OUTPUT'
      WRITE(99,*) 'CORE RECON POST_WDTH'  
      DO core=1,data_set%no_core_samples
        test_env_var=min_env_var
        post_width=0.0
        sum_pdf=0.0
        DO o=1,100
          post_width=post_width+data_set%core_rec_env_var_pdf(core,o)*&
            (test_env_var-data_set%core_rec_env_var(core))**2
          sum_pdf=sum_pdf+data_set%core_rec_env_var_pdf(core,o)
          test_env_var=test_env_var+env_var_interval
        ENDDO
        post_width=sqrt(post_width/sum_pdf)
        WRITE(99,990) core,data_set%core_rec_env_var(core),post_width,data_set%core_richness(core)
      ENDDO
  990 FORMAT(I5,2F10.3,I3)

      CLOSE(UNIT=99)

END SUBROUTINE DIAGNOSE_CORE
