function v = apply_quatrot(q,v)
v = [ zeros(size(v,1),1), v ];
v = multiply_quatrot( q, multiply_quatrot(v, conj_quatrot(q) ) );
v = v(:,2:end);
end
