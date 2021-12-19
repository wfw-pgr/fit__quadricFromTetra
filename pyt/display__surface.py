import numpy as np

# ========================================================= #
# ===  display__surface                                 === #
# ========================================================= #

def display__surface():

    x_, y_, z_ = 0, 1, 2
    
    # ------------------------------------------------- #
    # --- [1] make grid                             --- #
    # ------------------------------------------------- #
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ -1.0, 1.0, 21 ]
    x2MinMaxNum = [ -1.0, 1.0, 21 ]
    x3MinMaxNum = [  0.0, 0.0,  1 ]
    coord       = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "point" )
    import nkUtilities.load__pointFile as lpf
    coef    = lpf.load__pointFile( inpFile="dat/coef.dat" , returnType="point" )
    tetra   = lpf.load__pointFile( inpFile="dat/tetra.dat", returnType="point" )
    surf    = coef[0] * coord[:,x_]**2 + coef[1] * coord[:,y_]**2 \
        +     coef[2] * coord[:,x_]*coord[:,y_] \
        +     coef[3] * coord[:,x_]    + coef[4] * coord[:,y_] \
        +     coef[5]
    coord[:,z_] = surf
    
    import nkUtilities.save__pointFile as spf
    outFile   = "dat/surf.dat"
    spf.save__pointFile( outFile=outFile, Data=coord )
    
    # ------------------------------------------------- #
    # --- [2] convert surface into vtk              --- #
    # ------------------------------------------------- #
    import nkVTKRoutines.convert__vtkPolySurface as cvp
    outFile   = "dat/surf.vtp"
    cvp.convert__vtkPolySurface( Data=coord, outFile=outFile )

    # ------------------------------------------------- #
    # --- [3] convert point into vtk                --- #
    # ------------------------------------------------- #
    import nkVTKRoutines.scatter__vtkPoint as sct
    sct.scatter__vtkPoint( Data=tetra, vtkFile="dat/points.vtu" )
    
    
    return()


# ========================================================= #
# ===   Execution of Pragram                            === #
# ========================================================= #
if ( __name__=="__main__" ):
    display__surface()
