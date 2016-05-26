2/16/2010

Trifid
v1.0
Keith Wiley
kwiley@keithwiley.com
http://keithwiley.com

Note: throughout this document the Greek symbol 'mu' is often written with a lowercase Roman 'u'.

Note: throughout this document the term CCD is used to describe an electronic digital imaging sensor.  I believe that CMOS sensors are roughly equivalent in function within the level of description employeed here.

Note: the optical tools known to focus light are both the refractive lens and the concave mirror.  For simplicity, this document uses the word 'lens' throughout, but bear in mind that 'mirror' could be substituted virtually without error.

Sections:

    -- Introduction
    -- Sample Usage
    -- File Names
    -- Flags
    -- Noise Parameters File Description
    -- Warning Summary
    -- Denoising Flag Description
    -- Appendix 1: Scene Luminance Unit Conversion
    -- Appendix 2: Reflected lux to CCD incident PPFD Conversion
    -- Appendix 3: Stellar Focal Plane Illumination
    -- Appendix 4: The Absence of Lens Diameter from the Model
    -- References

--------------------------------------------------------------------------------
VERSION HISTORY
100222
v1.0 (not changing the name for this version)
-- Added mention of Weiner filter to README.

1002??
v1.0
-- First public release

####################################################################################################
####################################################################################################

INTRODUCTION
========================================

Trifid is a program which provides a noise model of CCD sensors and digital cameras.  The model more accurately represents the types of noise that are likely to occur in real world imaging than would generally be achieved by using a more naive approach, such as simply adding Gaussian noise to an image.  While the discussion of many of the parameters of this model focuses on CCD imaging, there is good reason to believe it will provide a comparable benefit for film noise since most sources of noise in either method of imaging map directly onto an analogous source of noise in the other method.  Although I have not objectively verified that this model produces more natural noise than a simpler approach, the results are nevertheless subjectively much more pleasing to my eye, and I am reasonably sensitive to the characteristics of a variety of natural types of image noise since I have spent many years capturing and processing images under extremely noisy conditions (long exposure astrophotography using low-end camera equipment).

One obvious direction for future work would be to more rigorously quantify the noise characteristics intrinsic to real imaging scenarios and to verify that this model really does produce results which are more accurate than a more naive approach.

The name Trifid comes from a popular astronomical object, M20, aka NGC 6514, aka the Trifid Nebula, which presented a fairly challenging subject back when I did a lot of astrophotography.  The results were always...well...noisy.  You can see my best attempts on my website.

--------------------------------------------------------------------------------
METHODS OF OPERATIONS

Trifid runs in a few different modes of operation, only one of which is the addition of noise to an image.  The full set of operating modes is as follows:

    -- Generate simulated noisy images.
    -- Denoise images using soft shrinkage and the steerable pyramid with four oriented filters.
    -- Generate image statistics between noisy or denoised images and an original baseline image.

See the -h flag for a full description of Trifid's modes of operation the flags and parameters it accepts as input.

--------------------------------------------------------------------------------
TERMINOLOGY

The following terms should be at least casually understood before attempting to comprehend this document or how Trifid works.

    -- Luminosity function: Description of the sensitivity of the human eye to different wavelengths.

    -- Luminous intensity: Power, weighted by the luminosity function, emitted by a light source in a particular direction per unit solid angle.

    -- Luminance: Luminous intensity per unit area in a given direction.  The total amount of light that passes through or is emitted from a particular area and falls within a given solid angle.  Often used to describe how bright a scene is.  It will be the driving parameter of overall scene brightness in this noise model.
           Common units:
               -- cd/m^2 (candela per square meter)
               -- lux, sometimes, not clear on how this represents luminance

    -- Luminous flux: Total luminous energy, weighted by the luminosity function, incident on or emitted by a surface.
           Common units:
               -- lm (lumen)

    -- Illuminance: Total luminous flux incident on or emitted by a surface per unit area.  Intensity of incident light, weighted by the luminosity function.
           Common units:
               -- lux (lm/m^2 lm lm is lumens)

    -- Photon Flux: The actual number of photons arriving within a given area over a given duration during a specific imaging session.
           Common units:
               -- photons

    -- Photon Flux Density: The flux rate per time.
           Common units:
               -- photons/m^2/s

    -- Photosynthetic Photon Flux Density: flux density confined to the frequency range [400nm, 700nm]
           Common units:
               -- ppfd (umols of photons/m^2/s
           Synonymous with PAR, Photosynthetic Aactive Radiation

--------------------------------------------------------------------------------
PROPER CONCEPTUALIZATION OF THE INPUT IMAGE

A new user might expect that the proper (and only) concept of Trifid's behavior is simply that it adds noise to an input image, i.e., it produces noisy output images from unnoisy (or less noisy) input images.  While this level of understanding will be sufficient for simple usage, a more rigorous understanding of how Trifid works will assist in advanced usage, namely, it will assist in adjusting the settings of the noise parameters description file.

In a more complex conceptualization, the scene's illuminance shining on the scene depicted in the input image is specified by a single parameter in the noise parameters file (the 'lux' parameter).  This parameter is generally assigned from a table of standard luminances of natural scenes (provided in this document).  For example, conventional outdoor daylight scenes are illuminated by about 10,000 to 25,000 lux.  The purpose of the input image is to indicate the spatial (pixel) and frequency (in three RGB bins) distribution of three inseparable properties of the scene within the simulated camera's field of view:

    -- surface-material reflectance, range: [0, 1], which linearly attenuates reflected luminance
    -- angle-of-reflection, range: [0, pi/2], which results in a cosine attentuation of reflected luminance
    -- shadow by surfaces residing between the original light source and the surface location in question

As a triplet, these three properties will be termed 'scene reflectance' through-out this document (not to be confused with the first property, material reflectance).  For a given pixel and RGB channel in the input image, the overall scene luminance and that pixel's scene reflectance will combine to indicate the subsequent luminance that the pixel will emit toward the camera, termed 'reflected luminance' (per-pixel-per-channel).  The reflected luminance will then pass through the lens optics to project a subsequent per-pixel-per-channel-luminance onto a photo-sensitive pixel of the CCD, termed 'sensor incident luminance'.

Please see the appendix for further elaboration of this discussion.

Note: We will assume that the angle-of-view represented in the scene perfectly fits the angle-of-view of the CCD and that the resolution of the CCD is identical to that of the input image.  Thus, the input image maps one-to-one onto the CCD and the resulting output image.

Note that while luminance and illuminance are generally expressed in units of lux, this model operates on a different representation, 'photosynthetic photon flux density', or ppfd.  This concept is explained in detail in the appendix.

--------------------------------------------------------------------------------
NOISE MODEL

Given a ppfd in photons/px/s (for a given pixel and channel), we must then convert to a digital measure of the corresponding pixel's brightness which will be represented in the output image.  At the same time a noise model will be applied to simulate the realistic sources of and variations in noise that can occur.

Note that while it is theoretically possible to process the ppfd-to-output-signal pipeline without actually adding any noise, the resulting output image would not be identical to the original input image.  The reason is two-fold.  First, as described in the appendix, the pixel values in the input image will have already been heavily processed to convert to ppfd incident at a CCD pixel.  Second, the ppfd will undergo yet further processing (even without the addition of noise) to generate a final output image.  Thus, the input and output images can never be assumed to be in perfect correspondence, even in the complete absence of the noise model.

We will need to know the quantization step size of the ADC in two places below, so let's calculate it first:

Ref: [2] I think.  Can't remember, can't find online preview of the relevant pages

    quantizationStepScaled = quantizationStep * fullWellCapacity
    where:
      quantizationStep = 1.0 / (adcRequiredBits^2 - 1)
      where:
        adcRequiredBits = the integer number of bits required to capture the CCD's dynamic range, i.e., the minimum integer satisying:
          6.02 * adcRequiredBits - 3.1876 = adcDynamicRange >= sensorDynamicRange
          where:
            sensorDynamicRange (dB) = 20.0 * log_10(sensorSNR)
            where:
              sensorSNR = fullWellCapacity / readNoise
              where:
                fullWellCapacity = an input parameter, the maximum integer number of electrons a single pixel of the CCD can store before saturation
                readNoise = an input parameter, the CCD's noise floor, the mean number of electrons generated per pixel in the absence of illumination and prior to gain-amplification

We will need to know the variance of signal-independent noise, i.e., variance of noise which will not vary by pixel or channel, i.e., variance which has no correlation with and is not a consequence of the scene.  Our initial version of this calculation (to be elaborated as we go) is the following:

Ref: [1: p. 270, Eq. 25]

    signalIndepVar = ampGainSquared * ampVarFactor + readNoiseVar + quantizationStepScaled^2 / 12
    where:
      ampGainSquared = ampGain^2
      where:
        ampGain = a unitless scalar, either an input parameter (user specified ISO) or is calculated automatically (auto-ISO)
      ampVarFactor = a unitless input parameter used to calibrate amplifier associated noise
      readNoiseVar = readNoise^2

However, quantization noise (that resulting from (quantizationStepScaled^2 / 12)) will occur implicitly when the analog-to-digital conversion occurs later.  Thus, we can simplify this calculation to the following:

Ref: [1: p. 270, Eq. 25]

    signalIndepVar = ampGainSquared * ampVarFactor + readNoiseVar

We can now investigate how to convert the ppfd (photons/px/s) for a particular pixel and channel to a digitized output value:

  -- Convert ppfd to idealElectronCount (loosely, one-to-one)
  -- Convert idealElectronCount to actual electron count:
       actualElectronCount = min(
         idealElectronCount * bias + meanDarkCurrentElectrons + hotPixel + meanAmpGlowElectrons,
         fullWellCapacity)
         where:
           bias = quantumEfficiency (per RGB channel) * pixelDependendentSensitivity, although pixelDependendentSensitivity is ignored in this model
           meanDarkCurrentElectrons = an input parameter, the mean number of electrons generated per pixel per second at the current temperature
           hotPixel = pixelDependentHotPixelSeverity * temperatureScalar * exposureTime * hotPixelFactor * fullWellCapacity
             where:
               pixelDependentHotPixelSeverity = a location-dependent severity in the range [0, 1]
               temperatureScalar = 1 at 20C and doubles every +6C (and halves every -6C)
               exposureTime = an input parameter in seconds
               hotPixelFactor = a unitless input parameter used to calibrate hot pixel associated noise
           meanAmpGlowElectrons = an input parameter, the mean number of electrons generated per pixel per second by the amplifier
  -- Calculate signalDepNoiseVar = ampGainSquared * ampVarFactor * actualElectronCount
  -- Calculate totalNoiseVar = signalDepVar + signalIndepVar
  -- Calculate noise = a random value from a normal distribution with mean 0 and variance totalNoiseVar
  -- Calculate finalElectronCount = idealElectronCount + noise, clipped to the range [0, fullWellCapacity]
  -- Digitize the electron count through the ADC:
       digitalVal = ((finalElectronCount + halfQuantizationStep) / quantizationStep) * quantizationStep
  -- Normalize and clip digitalVal:
       digitalValNormClip = digitalVal / fullWellCapacity, clipped to the range [0, 1]

The final value digitalValNormClip is stored in the output image.  Bear in the mind that if the input image represented a Bayer mosaic filtered view, then the resulting image will still represent a single-channel Bayer mosaic filtered view and must be processed through a Bayer demosaicing interpolation algorithm to generate a color image.

--------------------------------------------------------------------------------
PREEXISTING BAYER DEMOSAICING ARTIFACTS

Note that for the purpose of experimentation, the actual color input image most likely already suffered from Bayer mosaic related artifacting.  Namely, if the input image was originally captured with a digital camera it most likely represents the result of a Bayer demosaicing process since most color digital cameras use a Bayer or similar array to generate a color image.  Even if the image originated from film, it must have been digitized somehow, perhaps with a flatbed scanner, or alternatively, by simply having a digital photo taken of it.  In the latter case, the presence of a Bayer demosaicing step is obvious, but even in the case of a scanner-generated image, a CCD was most likely used to convert the print to a digital medium within the scanner, and the most common way to get color information into a CCD is once again with a Bayer or similar color array.  Consequently, it is dangerous to treat the input image as an uncorrupted description of luminosity.  Thus, the input image will have already been demosaiced once before.  The solution to this problem is simple however.  Bayer demosaicing is effectively an undersampling problem, in that the pixel resolution is half the color-sampling resolution (each of the three RGB colors being sampled at every other pixel location), the demosaicing process being nothing more than an interpolation to fill in the unsampled pixels in each of the three channels (btw, the literature on Bayer demosaicing algorithms is extensive and Trifid uses one the best, "alternating projections", although applies it with what I believe is a superior wavelet to that generally used in its implementation.).  Since the root of the demosaicing problem is undersampling, the solution is straight forward: downsample the input image.  One might assume that a factor of two is all that is sufficient since the original undersampling was by a factor of two.  However, since many demosaicing algorithms consider a wide neighborhood of pixels, there may be some further benefit (with quickly diminishing returns I imagine) in downsampling by more than two, perhaps four.  Note that the necessity of downsampling probably won't be a problem in an experimental setting since an input test image will have most likely been captured with a camera that produces on the order of megapixels and which will therefore be vastly reduced in size anyway, simply to speed up experimentation.  Cropping, of course, doesn't count since it does not downsample.

While I have made something of an issue of Bayer demosaicing artifact problem here, note that in my experience this is rarely much of a problem.

--------------------------------------------------------------------------------
PREEXISTING INPUT IMAGE NOISE

Related to the note on Bayer mosaicing above, also note that the input image has inherently already suffered from all the types of noise being subsequently modeled here.  This need not be a problem, especially for experimental applications, provided the experimenter keeps in mind that they are explicitly assigning ground truth to their input image, i.e., treating it as a perfect description of natural luminosity with zero noise or other type of error, regardless of the fact that this is clearly not the case.  The Bayer demosaicing problem was mentioned individually because it is the type of preexisting corruption most likely to produce noticeable and possibly confusing artifacting if the user does not bear it in mind carefully.

####################################################################################################
####################################################################################################

SAMPLE USAGE
========================================

The following examples show how to perform Trifid's basic operations.  In addition, please look at the QUICKSTART.sh shell script and the TUTORIAL.txt file for more information on using Trifid.

Simulate the CCD response to a scene described by a grayscale PGM image:
    % ./Trifid -N noiseParamsFile image.pgm

Simulate the CCD response to a scene described by a color PPM image:
    % ./Trifid -N noiseParamsFile image.ppm

Simulate the CCD response to a scene described by a color PPM image by converting the color scene to a grayscale Bayer matrix, applying the noise model to that image, and then demosaicing the Bayer matrix to produce the final color result.  This comes the closest to simulating a full digital camera pipeline:
    % ./Trifid -N noiseParamsFile -b image.ppm

Denoise an image using four-orientation steerable pyramid shrinkage:
    % ./Trifid -n -s 0A.01 image.ppm

Compare PSNR and Lab Delta Error statistics of a set of noisy images against a baseline (presumably unnoisy) image.
    % ./Trifid --rgbStat baselineImage.ppm noisyImage1.ppm noisyImage2.ppm noisyImage3.ppm

####################################################################################################
####################################################################################################

FILE NAMES
========================================

Output filenames are automatically generated from input filenames, suffixed with brief codes to indicate the full processing history in chronological order from left to right.  Codes corresponding to a single processing step are separated by an underscore '_' while parameters describing a particular step are separated by a hyphen '-'.  The following is a list of the possible codes and parameters:

    Codes        Parameters    Description
    -----        ------        -----------
    _b                         A color image was Bayered (to gray)
    _db                        A gray image was deBayered (to color)
                -nn                Nearest neighbor deBayering method
                -bl                Bilinear deBayering method
                -hue               Smooth hue deBayering method
                -ap                Alternating projections deBayering method
    _n                         The noise model was applied
                -<string>          Label specified in the noise parameters file
    _dn                        The image was denoised
                -<string>          The sequence of '-s' commands used to set thresholds
    _out                       Ubiquitous output suffix, appended to all files, stripped (with the extension) from all input files

Full example:
    RoswellEvidence_b_n-eveningHandheld_db-ap_dn-0A.02-1H1.2_out.ppm

Explanation of the example:
The original input image's filename was RoswellEvidence.ppm (substr "RoswellEvidence").  It was clearly a color ppm image because the first thing that was done to it was Bayerisation to a gray image (substr "_b").  Next followed the application of the noise model using a noise parameters file whose label parameter was "eveningHandheld" (substr "_n-eveningHandheld").  The resulting noisy gray image was then deBayered using alternating projections to produce a color image (substr "_db-ap").  At this point, the image must have been saved to a file named "RoswellEvidence_b_n-eveningHandheld_db-ap_out.ppm" because the run stream does not currently provide a means to pass an image from one command to the next.  On a subsequent run the resulting image was read in and denoised using two '-s' flags to indicate shrinkage thresholds (substr "_dn-0A.02-1H1.2").  The final result was saved to the filename shown in the example.

####################################################################################################
####################################################################################################

FLAGS
========================================

Trifid excepts the following flags:

    -h                  Display program help
    -v                  Verbose output
    -i                  Produce all intermediate images

    --testFT            Test Fourier Transform Output and Speed
    --testFTconv        Test Fourier Convolution
    --testColor         Test RGB XYZ LAB Color Routines

    -N                  Add simulated sensor noise
                            In:  1+ pgm/ppm
                            Out: 1+ pgm/ppm
    -b                  Convert RGB input to Bayer gray before processing
                            Only valid when following -N in the arg list
    -n                  Denoise images
                            In:  1+ pgm/ppm
                            Out: 1+ pgm/ppm
    -s [threshold]      Specify a shrinkage threshold (see README for syntax)
                            Only valid when following -n in the arg list
    --grayStat          Calculate grayscale image statistics
                            In:  1  pgm, original
                            In:  1+ pgm, result
    --rgbStat           Calculate RGB image statistics
                            In:  1  ppm, original
                            In:  1+ ppm, result

####################################################################################################
####################################################################################################

NOISE PARAMETERS FILE DESCRIPTION
========================================

The noise parameters file consists of several settings, one per row, in a required order.  Blank rows and rows designated as comments by the presence of a leading pound-sign ('#') will be ignored.  Furthermore, the only part of an uncommented row that counts is the first whitespace-delimited "word", usually a number, occasionally a string.  Everything else on the line is considered a comment and is ignored.  By convention, each row contains the following entries, all but the first being ignored comments:

    -- numerical parameter assignment
    -- the parameter's type (integer, floating point, string)
    -- the setting's legal range, square-brackets inclusive, parentheses exclusive
    -- a brief description of the setting

Keep in mind that everything after the first "word" (after the first space) is solely for user clarification and will be ignored by Trifid, but by convention should be left in place.

The following list describes the items in the noise parameters file in the order they must occur.

========================================
filenameLabel
--------------------------------------------------------------------------------
string    [length > 0]
--------------------------------------------------------------------------------
This parameter is a short string of your choosing (no spaces allowed) that will be appended to the output filename.  Assuming you specify unique labels in different noise parameter files, this feature will enable you to determine which noise parameters were applied to a particular image.

========================================
lux
--------------------------------------------------------------------------------
float    (0, inf]
--------------------------------------------------------------------------------
This parameter specifies the luminance (sometimes called illuminance) of the overall scene.  A table listing the luminosity of a variety of natural lighting conditions can be found in the appendix.

CCDs react to discrete units of electro-magnetic radiation (photons), specifically in the visible range, but scene luminance is indicated in lux.  Although the user does not have to understand the relationship between these two measures, such an understanding will assist in a more refined usage of this model.  Please see the appendix for a detailed description of this unit conversion.

========================================
minReflectance
--------------------------------------------------------------------------------
float    [0, 1)
--------------------------------------------------------------------------------
The scene reflectance corresponding to a pixel of the input image with minimum value, i.e., 0.  Minimum scene reflectance may very well not be 0 itself since virtually all material reflectances are greater than 0 (all materials reflect some light), virtually all surfaces have some nonzero angle-of-reflection, and virtually all surfaces receive some nonshadow-occluded light either via transparency, penumbra, or mult-bounce reflection.

Please read the section above on interpreting the input image.

========================================
maxReflectance
--------------------------------------------------------------------------------
float    (minReflectance, 1]
--------------------------------------------------------------------------------
The scene reflectance corresponding to a pixel of the input image with maximum value, i.e., 255 in an 8-bit image.  Maximum scene reflectance is rarely 1 since virtually no materal reflects light with 0 absorption and virtually no surface has zero angle-of-reflection (the light source would be coincident with the camera lens in such a case!).

Please read the section above on interpreting the input image.

========================================
temperatureC
--------------------------------------------------------------------------------
float    [-100, 50]
--------------------------------------------------------------------------------
"Room temperature" is, by definition, assigned precisely 20C in this model.  For every 6C increase in temperature there will be a two-fold increase in dark current and hot pixel intensity, and vs/va for a 6C decrease in temperature.  In real world situations, normal photography would probably occur between -10C and 30C.  Astronomical CCD cameras are frequently mechanically cooled to as low as -20C for low-end equipment and -100C for high-end scientific equipment.

========================================
pixelPitchMicrons
--------------------------------------------------------------------------------
float    (0, inf]
--------------------------------------------------------------------------------
This parameter sets the pitch of a pixel on the CCD, i.e., the length along one side.  This model assumes that pixels are square.  Although such an assumption is often not quite true, I chose to dispense with that particular headache in this model.  Note that this parameter is not the surface area, it is the linear dimension.  Don't get the two confused.  Assuming a pixel's entire surface area is available for photo-conversion (again, a simplification of actual CCDs which I feel is justified at this level of detail), the surface area of a pixel will then simply be the square of this parameter's value.  Typical values are 1.9u for high megapixel, small sensor cameras, common in midrange cameras (which pack a large number of pixels onto a small chip), 5.6u for first decade webcams (which used similar sized chips but with much fewer pixels), 8-9u for highend DLSRs, 9-11u for high-end astroimaging cameras, and 15-25u for observatory cameras.

If you have a specific camera in mind, you might be able to find the exact CCD model used by the camera, and subsequently find the exact pixel pitch.  However, this specific data is frequently not available.  Nevertheless one can often approximate the pixel pitch from the CCD's overall size and the pixel count.  Please see the appendix for a description of this calculation.

Please see the note in the fullWellCapacity section pertaining to well density and the associated warning which may be produced during a run.

Please see the note in the exposureSecs section pertaining to saturation and the associated warning which may be produced during a run.

========================================
fullWellCapacity
--------------------------------------------------------------------------------
int    [1, inf]
--------------------------------------------------------------------------------
There is an upper limit on the number of electrons that can be stored in a collection site (a pixel) of a CCD sensor.  This upper limit is literally a function of volume, namely the surface area of a pixel combined with the "well depth".  The effect of a larger full-well capacity is to increase the captureable dynamic range of the camera, and likewise to decrease the overall noisiness of the resulting images by permitting longer exposures relative to a given amplifier gain before saturation.  The same benefit is associated with increased pixel surface area of course (both pixel surface area and full-well capacity increase with the expense of a camera and both contribute to lower noise properties).  One serious source of "noise" which is not included in this model is "blooming" or "bleeding" in which saturated pixels literally overflow their electrons into adjacent pixels (like filling an ice-cube tray).  This is the reason why saturated points of light, i.e., stars, often appear as large dispersed disks or are associated with full columns of saturated pixels.  Typical values for this parameter are 30,000-50,000 for low-end webcams, 80,000-100,000 for conventional cameras, 100,000-120,000 for high-end cameras (DSLRs), and 750,000 for observatory cameras.  That said, some CCDs have very small pixels and consequently very small full well densities, e.g., 6000 for 1.86u pixels.

Note: Intuitively this would be a function of the volume of the pixel, but there is an upper limit on the possible volume due to photon absorption rates through silicon.  Consequently, surface area generally governs full well capacity, not volume.

Trifid will report the well density in e/u^2.  Pay attention to this value.  If it is outside a realistic value a warning will be shown.  You can avoid this warning by changing the fullWellCapacity or pixelPitchMicrons parameter.  Electrons per square micron tend to be 400-1800, lower for older or cheaper cameras, higher for newer or, more expensive cameras.

Please see the note in the readNoise section pertaining to adc bits and the associated warning which may be produced during a run.

Please see the note in the exposureSecs section pertaining to saturation and the associated warning which may be produced during a run.

========================================
readNoise
--------------------------------------------------------------------------------
float (0, fullWellCapacity)
--------------------------------------------------------------------------------
This parameter conglomerates all pre-gain-amplifier circuitry-associated noise in the CCD.  Such noise is often generated from multiple sources.  It is independent of scene luminance or exposure time (I'm a little unclear whether it is independent of temperature) and is generally indicated in a constant mean number of electrons generated per pixel. It determines the low-limit noise floor, below which no discernible signal can possibly be recorded.  Typical values tend to be in the range 10-100, although a readNoise outside this range is also possible for CCDs at the extremes (higher for older or cheaper CCDs, lower for newer or more expensive CCDs).

When Trifid is run, it will produce a fair amount of textual output in addition to generating the final noised image.  One section of the textual output will convey information about the CCD's signal-to-noise ratio and associated dynamic range which are calculated directly from the fullWellCapacity and the readNoise in the following way:

Ref: [3]

    ccdSnr  = fullWellCapacity / readNoise
    ccdDRdB = 20 * log_10(ccdSnr)

Trifid will automatically calculate the minimum number of bits that the analog-digital-converter must use to represent the CCD's dynamic range.  For a given number of bits, the ADC's dynamic range will be:

Ref: [2: p. 420, Eq. 11.20]

    adcDRdB = 6.02 * adcBits - 3.1876

Trifid will assign adcBits such that adcDRdB meets or exceeds ccdSnr, and the resulting adcDRdB will be reported back to the user during the run.  Pay attention to the value of adcDRdB as reported by Trifid, i.e., the ADC's required bits.  Typical values are 10-14 bits for midrange to highend cameras and up to 16 bits for astrophotography cameras.  If you are getting values below 10, you are not simulating a very realistic CCD (although 8-bit ADCs aren't unheard of).  Trifid will produce a warning if adcBits is outside the realistic range.  Correct this by increasing fullWellCapacity or decreasing readNoise.  If you are getting values above 12-14 for a midrange to expensive camera or above 16 for an astrophotography camera, you are not simulating a very realistic CCD.  Correct this by decreasing fullWellCapacity or increasing readNoise.

========================================
ampVarFactor
--------------------------------------------------------------------------------
float    [0, inf]
--------------------------------------------------------------------------------
This parameter is an arbitrary constant that works in concert with the ampGain parameter to indicate the amplifier noise severity.  Somewhere between .3 and .5 seems to be subjectively realistic in the current model.  

========================================
darkCurrentRoomTempElectronsPerPixPerSec
--------------------------------------------------------------------------------
float    [0, inf]
--------------------------------------------------------------------------------
This parameter is a direct measure of the dark current electron-generation rate, i.e., e/px/s, at room temperature (20C).  I an unsure whether this parameter should scale with pixel area.  If you don't want to use this feature in your noise modeling simply set this parameter to 0.

This parameter seems to vary greatly between CCD models.  Some resources have suggested 10,000-25,000 e/px/s at 20C but Sony's HAD (Hole-Accumulation-Diode) CCDs are notorious for very low dark currents and are popular CCDs, having been used in many noncheap (but still consumer-grade) webcams.  I would assume more recent and more expensive CCDs are even better.  Relatively expensive Sony HADs (considerably better than those used in webcams) have dark current rates of only 1-10 e/px/s at 20C.  I'm guessing the numbers of medium-range webcams are 10-100.  I have no idea about the CCDs used in digital cameras, but they're probably pretty good, also in the 1-100 range.  I'm just guessing though.

========================================
hotPixelFactor
--------------------------------------------------------------------------------
float    [0, inf]
--------------------------------------------------------------------------------
This parameter arbitrarily calibrates hot pixel intensity across the full spectrum of temperature and exposure settings.  A value of 1.0 seems to approximately match real world experience.  I suggest leaving it alone.  If you don't want to use this feature in your noise modeling either set this parameter to 0 or set the hotPixelMask parameter to "NONE".  Either option will disable hot pixel noise.

Including hot pixel noise in the model requires a hot pixel mask image in .pgm format with the same dimensions as the input image.  Obviously, the best results will be achieved if the mask accurately represents the hot pixel distribution of a real CCD.  The 512x512 image I have provided (hotPixelMask512.pgm) seems to be a subjectively satisfying approximation, but it was not generated in a rigorous manner.

========================================
hotPixelMask
--------------------------------------------------------------------------------
string    ["NONE" or pgm_file_name]
--------------------------------------------------------------------------------
The name of a pgm file that represents the hot pixel pattern for a particular CCD.  A subjectively realistic example is provided with this project, named 'hotPixelMask512.pgm'.  Note that this image must have the same dimensions as the input image.

Please see the note in the hotPixelFactor section about this file.

========================================
quantumEfficiencyRed
quantumEfficiencyGreen
quantumEfficiencyBlue
--------------------------------------------------------------------------------
float    [0, 1]
--------------------------------------------------------------------------------
The quantum efficiency of the sensor would generally be represented as a smooth function across the visible spectrum.  However, in this model it is simplified to an analogy to the Bayer mosaic channels in the pseudo-wavelength bands of 400-500nm, 500-600nm, and 600-700nm.  The corresponding color channels will be reduced in efficacy by these QE parameters.

Note: Without the '-b' flag, these parameters have no effect on the simulation, even when applied to color images.

The quantum efficiencies of real CCDs vary greatly.  Some are weaker in the blue and stronger in the red, others vs/va.  Some, although not many, are relatively flat across the visible spectrum.  Examples of realistic sets include (.4R, .6G, .7B) or (.7R, .6G, .4B).  The consequence of varying QE by color -- that is, by wavelength -- is first, a color shift in the resulting image, and second, more noise in the lower QE colors after the color shift is corrected through a scaling that inverts the QE.  Such a scaling will enhance the lower QE colors to correct the color shift, thus amplifying their respective noise.

Note that if all three quantum efficiencies are equal then this parameter will have been effectively removed from the model since all pixel values will simply be scaled by the same amount.  Consequently, if you set the three quantum efficiencies to the same value, you may as well set them all to 1.0, which numerically removes them from the simulation.  In fact, ideally, at least one quantum efficiency could be 1.0 in all simulations even though this is not representative of true CCDs, given their linearity.  If you wish to remove this behavior from the model, simply set all three parameters to 1.0.

========================================
ampGlowElectronsPerPixPerSec
--------------------------------------------------------------------------------
float    [0, inf]
--------------------------------------------------------------------------------
The amplifier on a CCD chip produces actual photons which can then be detected by the CCD.  The amplifier is usually located just off one corner of the sensor.  If this feature were properly implemented, there would be a 1/radius dropoff in amp glow intensity w.r.t. one corner of the chip, but in the current model such variation is not included.  Amp glow is not generally a problem because, in the case of conventional photography, exposures are kept brief enough (fractions of a second) that amp glow does not accumulate an appreciable amount, and in the case of professional long-exposure photography (like astroimaging) the amplifier is powered down during the exposure (it is only needed during readout).  Amp glow is only a problem in the case of long exposures taken with cameras that do not power down the amplifier during such exposures.  This can occur when using long-exposure modified webcams, although slightly more advanced modifications permit the same amplifier power-down associated with high-end cameras.  Since this feature is not fully implemented (w.r.t. to the 1/radius dropoff) and since amp glow is not a noticeable problem either in short duration conventional photography or in long exposure photography with amp-power-down cameras, I recommend leaving this parameter set to 0 and effectively removing it from the model.

========================================
exposureSecs
--------------------------------------------------------------------------------
float    [0, inf]
--------------------------------------------------------------------------------
This parameter is typically in the range .0005s-.05s for daytime photography with a normal camera with no tripod, or for astroimaging of bright planets or the moon with a telescope.  It is typically in the range .05s-10s for normal photography with a tripod.  It is typically in the range 1s-100s for amateur astrophotography of deep sky objects (nebulae and galaxies) with an amateur telescope.  It is typically in the range 100s-1000s for professional astrophotography of deep sky objects.

Trifid will produce a warning if the exposure time exceeds the theoretical saturation time before the addition of noise.  If you get this warning, you need to decrease the photon flux or increase the saturation point.  This affect can be achieved in many ways: by decreasing lux, decreasing maxReflectance, increasing pixelPitchMicrons, increasing fullWellCapacity, decreasing exposureSecs, or increasing focalRatio.  However, a photographer would generally only have control over the focalRatio and exposureSecs parameters (within the ranges provided by their equipment).  An astrophotographer would might only have control over the exposureSecs depending on their equipment.

========================================
ampGain
--------------------------------------------------------------------------------
float    [-1][1, inf]
--------------------------------------------------------------------------------
Setting this parameter to -1 is a flag to indicate that auto-gain should be performed such that good saturation is achieved.  Doing so is analogous to auto-ISO on a digital camera in which the ISO setting directly controls the CCD amplifier gain.  Note that this parameter is specified as a literal scalar, not an ISO value as would conventionally be used when adjusting the corresponding setting on a camera.  This is because assigning an ampGain-to-ISO mapping to a particular camera is a rather arduous process involving real-world tests and comparison with film-established baselines.  Furthermore, the slightest change to any parameter of would require a new ampGain-to-ISO mapping to be calculated.

========================================
focalRatio
--------------------------------------------------------------------------------
float    (0, inf]
--------------------------------------------------------------------------------
The terms "focal ratio", "lens number", "f-number", and "f-stop" are all synonymous and all refer to the ratio of a lens's focal length to its diameter.  Thus follows a well-established relationship:

Ref: [11]

    N = f / D

where N is the focal ratio, f is the focal length, and D is the working aperture diameter.  Obviously, given any two of these parameters, the third is uniquely determined.  In conventional photography with a nonzooming lens, the photographer only has explicit control over the f-stop, i.e., the focal ratio.  While the lens's focal length never changes (it is a physical property of the curvature of the glass) the diameter does effectively change in that the diaphragm reduces the lens's working diameter from its full size to some fraction thereof.  The size of the diaphragm opening is controlled by the f-stop setting.  In a slightly more complicated scenario, the use of a zoom lens permits the photographer to control the parameter f, the focal length, as well, by modifying the relative distances of the various lens components.  Things are quite different in astrophotography however.  Generally all three parameters are fixed.  D is simply the width of the primary lens (or more frequently a mirror), f is a property of the curvature of the lens and N follows implicitly from f / D.  Sometimes astrophotographys introduce additional optics into the light path such that the overal focal ratio is altered.  Eyepieces and focal reducers are common examples of such additional pieces of equipment.

The effect of this parameter is that the illumination at the focal plane will reduce by the square of the focal ratio.  Thus, if you double the focal ratio between two experiments, the illumination will be cut by a factor of four.  Other parameters will have to be adjusted to compensate, conventionally exposure time and amp gain.

Please see the appendix for an important note on how N, f, and D are used by Trifid.

####################################################################################################
####################################################################################################

WARNING SUMMARY
========================================

This section lists the warnings that may be generated by Trifid and prescribes possible solutions, along with a few notes.

--------------------------------------------------------------------------------
ADC required bits (Xb) is unrealistically low (below 10b).

Solutions:
    -- Increase fullWellCapacity
    -- Decrease readNoise

--------------------------------------------------------------------------------
ADC required bits (Xb) is unrealistically high (above 16b).

Solutions:
    -- Decrease fullWellCapacity
    -- Increase readNoise

--------------------------------------------------------------------------------
Well density (Xe/u^2) is unrealistically low (below 400e/u^2).

Solutions:
    -- Increase fullWellCapacity
    -- Decrease pixelPitchMicrons

--------------------------------------------------------------------------------
Well density (Xe/u^2) is unrealistically high (above 1800e/u^2).

Solutions:
    -- Decrease fullWellCapacity
    -- Increase pixelPitchMicrons

--------------------------------------------------------------------------------
Exposure time (Xs) exceeds no-noise saturation time (Xs).

Solutions:
    -- Decrease lux
    -- Decrease maxReflectance
    -- Increase pixelPitchMicrons
    -- Increase fullWellCapacity
    -- Decrease exposureSecs
    -- Increase focalRatio

####################################################################################################
####################################################################################################

DENOISING FLAG DESCRIPTION
========================================

When run in denoising mode (use the -n flag to engage this mode), the -s flag is used to assign shrinkage thresholds.  The interface provides full control over the individual thresholds w.r.t. color and pyramid level.  However, thresholds cannot be indicated on a per-orientation basis.  The general syntax for a threshold parameter is one of:

    -s [color][A|H][threshold]
    ...or...
    -s [color][B][level][threshold]

where...

    [color] is one of:
        '0'          Set all three red, green, and blue channels
        '1'          Red
        '2'          Green
        '3'          Blue

    [A|H][B] is one of:
        'A'          Set all thresholds at all levels for the indicated color
        'H'          Set the H0 level threshold
        'B'          Set the B threshold at some level

    [level] is one of:
        '0'-'9'      Set the B threshold at the indicated level

    [threshold] is:
        float        Set the threshold value, decimal point or leading/trailing zeros optional

Crucially, notice that the syntax differs for B relative to A and H in that B requires that a level parameter be specified and A and H require that no level parameter be specified.  If used incorrectly, the presence or absence of the level parameter will become confounded with the threshold parameter, so be careful.

If multiple -s entries are provided on the command line, those appearing later in the argument list override those appearing earlier in the list, so specify thresholds from general to specific, as demonstrated below.

Examples:

Set all thresholds for the entire pyramid to .1:
    -s 0A.1

Set all red thresholds to .1:
    -s 1A.1

Set all H thresholds to .1:
    -s 0H.1

Set green B level 2 threshold to .1:
    -s 2B2.1

Set all thresholds to .1, then override all B level 1 thresholds to 5.0:
    -s 0A.1 -s 0B15

Set all thresholds to .1, then override all B level 1 thresholds to 5.0, then override red B level 1 threshold to 3.4:
    -s 0A.1 -s 0B15 -s 1B13.4

####################################################################################################
####################################################################################################

APPENDIX 1: SCENE LUMINANCE UNIT CONVERSION
========================================

This section elaborates on the topics of scene luminance (lux) and pixel illuminance (photon flux density).

The method of describing luminance which is most applicable to a study of CCD imaging is Photosynthetic Photon Flux (PPF), i.e., the discrete number of photosynthetically active photons (400nm-700nm) recorded by a given pixel of the CCD, or photons/area (photons/pixel).  Related to PPF, there is the similar measure Photosynthetic Photon Flux Density (PPFD), also called Photosynthetically Active Radiation (PAR) which indicates not the number of photons recorded, but rather the photon flux *rate*, i.e. photons/area/time (photons/pixel/second).  Although the names of these metrics clearly suggest an intended application to botany, they are also exactly what we are interested in when describing image capture for the purpose of human vision since the range of photosynthetic wavelengths almost perfectly corresponds to the range of human-visible wavelengths.  We use these terms and their common units of measurement because there is a substantial quantity of accessible data available on the associated botanical topics which can be applid to the design of imaging scenarios (noise parameter files).

PPFD is generally indicated in units of discrete photons as umols/m^2/s (also called microeinsteins/m^2/s) (remember that a lower case 'u' is an ASCII convention for a "mu" symbol, the standard prefix for "micro").  Recalling that one mole is Avogadro's number (6.0221415e23) we conclude that PPFD, as indicated in umols/m^2/s corresponds to 6.0221415e17 photons/m^2/s.

While PPFD is the unit of measure we can most easily process at the CCD, descriptions of natural scene luminance are generally provided in units of lux instead.  The following table lists the luminance for natural scenes in lux:

Ref: [4,5,6]

    -- Sun overhead                      130 klux, 50 klux summer Britain
    -- Full daylight (not direct sun)    10 klux - 25 klux, perhaps 50 klux
    -- Overcast day                      1 klux, 45 klux summer, 20 klux winter, 5 klux Britain
    -- Very dark overcast day            100 lux
    -- Sunrise/Sunset                    10 lux
    -- Twilight                          1 - 10 lux
    -- Deep twilight                     1 lux
    -- Full Moon overhead                .3 lux (direct, I believe, not scene reflectance)
    -- Cloudy moon:                      .1 lux
    -- Young/quarter moon:               .01 lux
    -- Total starlight + airglow         .002 lux
    -- Total starlight, overcast         .0001 lux

Note that these are measures of the light source, not of the objects being view in the scene.  Light received by the camera will potentially be significantly attenuated by the reflectance of materials in the scene, thus the following table of material reflectances:

Ref: [13,14]

    -- Snow                             90%
    -- Wall (painted white)             85%
    -- Sand                             50%
    -- Parking lot with cars            40%
    -- Concrete                         35%
    -- Grass                            30%
    -- Brick                            25%
    -- Black(asphalt,empty parking lot)  5%

Furthermore, various surfaces in the scene will reflect the original light source toward the camera with some angle of reflection such that steeper angles will attenuate the light according to a cosine function ("glancing" surfaces will not reflect light toward the camera as effectively as "face-on" surfaces).  Finally, various locations in the scene (pixels) may reside in the occluded shadow (umbra or penumbra) or transparent shadow of other surfaces.  Thus, the reflectances indicated in the table above represent an upper-bound for the materials in question.  The values stored in the input image will, by definition, represent an inseparable blend of these three properties (material reflectance, angle-of-reflection, and shadow), termed 'scene reflectance' throughout this document (on a per-pixel-per-channel basis).

Following a reduction in illuminance according to scene reflectance, the conversion from lux to PPFD is a straight-forward multiplication by .0185.  Thus, the following table represents the luminance for natural scenes (as listed above) in PPFD:

Ref: [7,8,9]

    -- Direct sun overhead               2405 umol/m^2/s, 925 summer Britain
    -- Full daylight (off-angle)         185 - 463, 925
    -- Overcast day                      18.5, 833 summer, 370 winter, 92.5 Britain
    -- Very dark overcast day            1.85
    -- Sunrise/Sunset                    .185
    -- Twilight                          .0185 - .185
    -- Deep twilight                     .0185
    -- Full Moon overhead                .005
    -- Cloudy moon:                      .002
    -- Young/quarter moon:               .0002
    -- Total starlight + airglow         .00004
    -- Total starlight, overcast         .000002

PPFD, combined with per-pixel-per-channel scene reflectance, and then processed through the optics and exposed for a prescribed duration, will ultimately determine the number of photons recorded by the pixel in question, i.e., the photon flux at the pixel.  This photon flux will then be subject to the noise model prior to the assignment of a final output value of the pixel.

####################################################################################################
####################################################################################################

APPENDIX 2: REFLECTED LUX TO CCD INCIDENT PPFD CONVERSION
========================================

One of the most important functions in the model is the conversion of reflected luminance in lux to a measure of CCD incident luminance as ppfd (ppfd isn't a unit per se, but it is always expressed in umol/m^2/s).  The conversion is peformed in the following manner:

For a given pixel...
    -- Start with scene lux.
    -- Reduce scene lux to lens transmission (account for lens reflectance and absorption).
    -- Reduce lens-exit lux by the square of the lens focal ratio (spread the light out at the image plane).
    -- Convert luminance in lux (lumens/m^2) to luminance in CCD pixels (lumens/CCD-pixel).
    -- Convert luminance from lumens to photons (to ppfd, i.e., photons/px/s)

The result of the pipeline described above is the number photons incident at a CCD pixel for a perfectly reflecting material.  Three steps remain prior to application of the noise model.  The first is to discard two of the three RGB channels for the pixel in question.  This step is implicit in the conversion of the input image via a Bayer mosaic operation.  The second step is to reduce the ppfd by the scene reflectance indicated in the input image (per-pixel-per-channel).  The third and final step is to expose the pixel to the final ppfd for the proper exposure duration.  The result will be the number of photons incident at the CCD pixel in question.  This number of photons will represent the pixel and channel dependent input to the noise model.

####################################################################################################
####################################################################################################

APPENDIX 3: STELLAR FOCAL PLANE ILLUMINATION
========================================

Please recall the note from the top of this document pertaining to lens/mirror equivalency.

Unlike the objects imaged in normal photography, in astrophotography there are two major classes of objects, stars and everything else.  What makes stars so special is that they are effectively point sources of light, they have no discernable width or area.  All other astronomical objects (planets, moons, asteroids, comets, nebulae, galaxies) have a measurable width and area, which is to say that they have a 'surface brightness' as astronomers call it.  This concept is meaningless for stars.  Consequently, the illuminance produced by an astronomical object varies with different properties of the optical system depending on whether that object is a star or something else.  In the case of all objects except stars, optical equipment works as we would expect, i.e., as it works in conventional photography, in which the illuminance imparted by an object at the focal plane is a function of the focal ratio of the imaging system.  In the case of stars, however, the illuminance produced at the focal plane is a function not of the focal ratio, but solely of the primary objective diameter (or area actually).  This distinction is clearly not made by this model since the model has no way of knowing which pixels in the input image correspond to stars.  This problem represents an area of potential extension or improvement to the model in the future.  An additional parameter would be required, lens diameter.  This point is expanded on further in the next section.

####################################################################################################
####################################################################################################

APPENDIX 4: THE ABSENCE OF LENS DIAMETER FROM THE MODEL
========================================

Please recall the note from the top of this document pertaining to lens/mirror equivalency.

Astrophotographers will be quick to notice the absence of an important parameter of an imaging system from this model, the diameter of the optics.  In the case of astrophotography, the diameter of interest is that of the primary objective, usually a large concave mirror.  Diameter is of critical importance in astrophotography because, when imaging objects of extremely limited illuminance, it is helpful to gather as many photons as possible.  Astrophotographyers jokingly refer to large telescopes as "light-buckets" for this reason.  Therefore, the question arises: how can this model justifiably discard lens diameter as an meaningful parameter?

The important thing to remember is that regardless of diameter, focal plane illuminance remains solely a property of the focal ratio (excepting stars for the current discussion, see the previous Appendix topic for more information about stars).  In other words, given two imaging systems with different lens diameters but equivalent focal ratios, both systems will produce the same illuminance at the focal plane, which is to say that a given exposure time will yield the same brightness of captured image in both systems.  How can this be so when we know for a fact that larger telescopes capture more light?  The key is to realize that in the example just given, the larger of the two imaging systems produces a larger image at the focal plane.  Think about it.  For the larger system to have the same focal ratio as the smaller system, it must have a longer focal length.  The focal plane is consequently further away from the lens, so the projected image is larger...but it is also dimmer because the light has been spread out...but it is also brighter because the diameter is larger and captures more light.  The wider diameter and longer focal length cancel one another out perfectly such that the focal plane illuminance is unchanged between the two imaging systems in question.

When we say that a larger telescope produces brighter images, what we really mean is that it produces brighter images at a given projected image size.  What would happen if we took the example from above and altered it so that the projected images were the same size?  First, how is this achieved?  We do this by decreasing the focal ratio of the larger telescope, i.e., decreasing the focal length.  The focal plane is pulled closer, the light is less spread out, and so the image is smaller.  Naturally, the image will now be brighter as well.  Thus, for an image of a given size (achieved with different focal ratios) the larger telescope produces a brighter image.  Bear in mind that most telescopes don't actually permit adjustment of the focal length (without additional optics, e.g., eyepieces and focal reducers); this is merely a thought experiment to make a point.

We can now ask how (why) lens diameter can be discarded from the model.  The model includes a focal ratio setting, but no diameter or focal length setting.  Recall the following equation from earlier in the document describing the relationship between focal ratio, focal length, and diameter:

Ref: [11]

    N = f / D

where N is the focal ratio, f is the focal length, and D is the working aperture diameter.  Given any two of these parameters, the third is uniquely determined.  However, the equation has two unknowns.  Only the focal ratio is specified as a parameter of the model.  The focal length and diameter are both unspecified.  By implication, we cannot meaningfully describe the size of the projected image.  But there is another parameter of the model, pixelPitchMicrons, which states the size of a single pixel of the CCD.  Given the dimensions of the input image and the implicit assumption that the input image specifically defines the angle-of-view of the scene (and the projected image), we can actually calculate the size of the image at the focal plane.  It is simply pixelPitchMicrons multiplied by the input image width (or height).  So on the one hand the image size is unconstrained, but on the other hand it is explicitly calculable.

The only logical interpretation of these observations is that both the focal length and the diameter are automatically scaled internally by the model so as to conform to the constraints specified by the parameterized focalRatio, pixelPitchMicrons, and input image resolution.  If diameter is effectively assigned automatically for us, we don't need it as a parameter and can ignore it.  Note that if you look at the code you will not actually find the lens diameter and focal length being explicitly calculated since they are unnecessary for the internal calculations of the model.  This touches on the previous Appendix topic which discusses how stars ought to be properly handled by the model.  Stars are not properly handled, but if they were, lens diameter would be required as a parameter.

####################################################################################################
####################################################################################################

REFERENCES
========================================

[1] (CCD noise model, referenced copiously in the code)
Radiometric CCD Camera Calibration and Noise Estimation,
Glenn E. Healey and Raghava Kondepudy,
IEEE Transactions on Pattern Analysis and Machine Intelligence,
Vol 16., No. 3, pp. 267-276, Mar 1994.
http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=276126

[2] (CCD noise model, referenced copiously in the code)
Analog & Digital Signal Processing,
H. Baher,
John Wiley & Sons, 2nd ed., 2001.
http://books.google.com/books?id=PqT9LHo7NO4C&pg=PP1&dq=baher+Analog+%26+Digital+Signal+Processing

[3] (CCD noise description)
http://www.qsimaging.com/ccd_noise.html

[4] (luminance of natural scenes and astronomical objects)
http://stjarnhimlen.se/comp/radfaq.html#10

[5] (luminance of seasons)
http://www.sunpipe.co.uk/technical/lux.php (out-of-date, seems to have switched to monodraught, listed below)
http://www.monodraught.com/technical/lux.php

[6] (luminance of natural scenes)
http://www.use-ip.co.uk/datasheets/lux_light_level_chart.pdf

[7] (PAR and conversion to/from lux)
http://openwetware.org/images/e/e8/Conversion_lux.pdf

[8] (minimal description luminance of natural scenes, also shows lux/ppf conversion)
http://www.apogee-inst.com/conv_lux.htm

[9] (lux/ppf conversion)
http://www.allcat.biz/mesurez/anglais/default/news.php

[10] (camera optics)
http://en.wikipedia.org/wiki/Angle_of_view

[11] (camera optics)
http://en.wikipedia.org/wiki/F-number

[12] (excellent explanation relating object brightness to image brightness through an optical system which is critical for calculating image plane illuminance after lens transmission.  Requires lens magnification as an input parameter.)
http://people.rit.edu/andpph/text-illuminance.html

[13] (camera optics)
http://en.wikipedia.org/wiki/Magnification

[14] (luminance of natural scenes, surface reflectance of common materials, description of camera optics, e.g. f-stop)
http://www.scansourcesecurity.com/MicroSites/ScanSourceSecurity/ipcenter/files/Cameras.pdf

[14] (surface brightness)
http://en.wikipedia.org/wiki/Albedo
