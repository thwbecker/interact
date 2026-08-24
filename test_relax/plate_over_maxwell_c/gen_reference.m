% gen_reference.m
%
% generate reference tables by running the original MATLAB codes
% (under Octave, together with the local roots.m wrapper that
% restores MATLAB's root ordering) at a set of test stations;
% used by run_verification.sh to check the C port
%
% each output file has one row per station: x y east north up

x=[-150 -50 0 50 150 -100 100 30 -30 250];
y=[ -80  40 60 10 -60  120 -120 200 -200 0];
xloc=[x;y];
H1=20; H2=40; nu=0.25;

function writeref(fn,x,y,U)
  f=fopen(fn,'w');
  fprintf(f,'%g %g %19.12e %19.12e %19.12e\n',[x;y;real(U)]);
  fclose(f);
end

% layer over halfspace: strike-slip, displacements and velocities
m=[200 20 20 90 0 0 0 1 0 0];
U = Plate_over_Maxwell_Layer_over_Halfspace_Displacements(m,xloc,H1,H2,nu,5,25,100);
writeref('ref_L_disp_ss.txt',x,y,U);
U = Plate_over_Maxwell_Layer_over_Halfspace_Velocities(m,xloc,H1,H2,nu,5,25,100);
writeref('ref_L_vel_ss.txt',x,y,U);

% dip-slip on a 30 degree dipping fault
m=[200 20 20 30 0 0 0 0 1 0];
U = Plate_over_Maxwell_Layer_over_Halfspace_Displacements(m,xloc,H1,H2,nu,5,25,100);
writeref('ref_L_disp_ds30.txt',x,y,U);

% tensile, oblique strike, offset fault position
m=[100 30 25 60 40 15 -10 0 0 1];
U = Plate_over_Maxwell_Layer_over_Halfspace_Displacements(m,xloc,H1,H2,nu,10,25,100);
writeref('ref_L_disp_ten.txt',x,y,U);

% mixed slip with equal relaxation times (ill conditioned; see README)
m=[200 20 20 70 20 5 5 1 0.5 0];
U = Plate_over_Maxwell_Layer_over_Halfspace_Displacements(m,xloc,H1,H2,nu,5,25,25);
writeref('ref_L_disp_mix_eq.txt',x,y,U);

% plate over Maxwell halfspace, mixed slip, three times
m=[200 20 20 90 0 0 0 1 0.3 0.1];
[Ue,Un,Uv]=Plate_over_Maxwell(m,xloc,20,1,[1 5 10],25);
for it=1:3
  U=[Ue(:,it)';Un(:,it)';Uv(:,it)'];
  writeref(sprintf('ref_M_t%d.txt',it),x,y,U);
end
