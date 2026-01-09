function[X,W] = triangulate_simple_og(x,L,P,TH)

[Nt,Np,Nc] = size(x);
Np = Np/3;
X = NaN*zeros(Np,3,Nt);
W = zeros(Np,Nt);

for t = 1:Nt
    for p = 1:Np
        cam_ind = [];
        temp = [];
        for c = 1:Nc
            if L(t,p,c)>TH
               cam_ind = [cam_ind c];
               temp = [temp [x(t,3*(p-1)+2,c); x(t,3*(p-1)+1,c); 1]];
            end
        end
        if numel(cam_ind)>1
           temp3D = ls_triangulate(temp,P(cam_ind)); 
           X(p,:,t) = temp3D(1:3)';
           W(p,t) = numel(cam_ind);
        end
    end
    if mod(t,100) == 0
       disp(sprintf('Triangulated %s frames out of %s', num2str(t),num2str(Nt)));
    end
end