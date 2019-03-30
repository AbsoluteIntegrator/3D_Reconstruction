
figure(1)
subplot(1,2,1)
imagesc(population(:,1:sequence_length))
title('Height Offset Map')
xlabel('Columns')
ylabel('Population')
colorbar

subplot(1,2,2)
imagesc(population(:,(sequence_length+1):(2*sequence_length)))
title('Atom Difference Map')
xlabel('Columns')
ylabel('Population')
colorbar

figure(2)
clf(2)

subplot(5,1,1)
plot(best_mass_log)
grid on
xlabel('Generation')
ylabel('Atoms')
title('Particle Mass of Best Active Model')

subplot(5,1,2)
plot(energy_per_atom_log,'r-')
grid on
xlabel('Generation')
ylabel('eV/atom')
title('Energy per Atom of Best Active Model')

subplot(5,1,3)
plot(model_probability_log)
grid on
xlabel('Generation')
ylabel('Probability')
title('Model Probability of Best Active Model')

subplot(5,1,4)
plot(cost_function_log,'g-')
grid on
xlabel('Generation')
title('Cost Function of Best Active Model')

subplot(5,1,5)
plot(inbreeding_log)
hold on
plot(child_diversity)
hold off
grid on
ylim([0 1.05])
xlabel('Generation')
ylabel('Uniquenes')
title('Model Uniqueness')

pause(0.1)
