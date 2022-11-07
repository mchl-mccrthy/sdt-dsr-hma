% Pixel shading

% This function was originally written by Neil Arnold with the name 
% f_whether_shaded, then modified very slightly by Michael McCarthy, to 
% determine shading when the shading of surrounding pixels 
% is unknown.

function [shaded] = f_whether_shaded_pixel(Cellelev,Pnt0,~,...
    DistStepVec,DistStep,ScaleZ,TanSolAlt,Dimx,Dimy)

% f_whether_shaded To find the shading for point 'pnt' in direction 'Vec'
%
%     Returns: Shade (0=shade, 1=sun)
%     By NSA, CUGD, April 1994
%     Matlab port 2014

%     Parameters and definitions
%    ----------------------------
%     Parameter      (Maxx=201,Maxy=301)
%
%     Passed Parameters
%    -------------------
%     Dimx                ! DTM X dimension             (Val)
%     Dimy                ! DTM Y dimension             (Val)
%     Pnt0                ! Start point for transect    (Val)
%     Vec                 ! Direction                   (Val)
%     DistStepVec         ! Step vector                 (Val)
%     DistStep            ! Sample rate from transect   (Val)
%     ScaleZ              ! Scaling for Dtm Heights     (Val)
%     TanSolAlt           ! Tangent of solar altitude   (Val)
%     Cellelev(Maxx,Maxy) ! DTM                         (Val)
%     Shadeimg(Maxx,Maxy) ! Shaded image            (Val/Res)
%
%     Get Height Of Initial Point
      PixCl  = fix( real(Pnt0) );
      PixRw  = fix( imag(Pnt0) );
      Pnt0Hgt = ScaleZ * Cellelev(PixCl,PixRw);
%
%     Set Initial Position
      Pnt        = Pnt0;
      Dist       = 0;
      Done       = false;
      Shade      = 1; % Pixel is initially sunny
%
%     For Each Point On Transect
      while (~Done)
         Pnt    = Pnt + DistStepVec;
         Dist   = Dist + DistStep;
%
%        Get DTM Pixel To Sample (Use Nearest Neighbour)
         PixCl = fix( real(Pnt) );
         PixRw = fix( imag(Pnt) );
%
%        Check For Off Edge Of DTM
         Done=((PixCl <= 0   ) | ...
     +         (PixRw <= 0   ) | ...
     +         (PixCl >= Dimx) | ...
     +         (PixRw >= Dimy));
%
%        If Not Off Edge
         if (~Done)
%
%           Compute Tangent Of Shading Angle For This Point
            Hgt     = ScaleZ * Cellelev(PixCl,PixRw);
            SolHgt  = Dist * TanSolAlt;
%
%           If new point is higher than height of sun, pnt0 is shaded...
            if (Hgt-Pnt0Hgt > SolHgt)
               Shade=0;
               Done=true;
            end
%
         end
%
%     If new point is not higher, pixel remains sunny and keep on walking
      end
%
%     Return
%     --------
      shaded = Shade;
end