
% Plot the RDF:
subplot(2,2,1)
hist(angstrom_distances(angstrom_distances>0),512)
xlabel('Angstrom')
title('RDF')

subplot(2,2,2)
plot(host_electron_density)
xlabel('Site Number')
title('Host Electron Density')

subplot(2,2,3)
plot(embedding_function)
xlabel('Site Number')
title('Embedding Function')

subplot(2,2,4)
plot(pair_potential_term)
xlabel('Site Number')
title('Pair Potential Term')
