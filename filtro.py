" Rotina teste para ilustrar a filtragem de uma série temporal na banda entre 3 e 30 dias"
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

''''
  TODO: FAZER REAMOSTRAGEM PRA 30 DIAS, PRA MELHORAR O DESEMPENHO NOS EXTREMOS
'''

def filtra_dados(nivel_mar, tempo, nome_serie, fig_folder, metodo, modelo = False):
  if metodo not in ['band', 'low', 'high']:
    metodo = input('ESCREVA low, high, ou band para selecionar o metodo.')
  elif metodo ==  'band':
    b, a = passa_banda(fig_folder, modelo)
  elif metodo == 'low':
    b, a = passa_baixa(fig_folder, modelo)
  elif metodo == 'high':
    b, a = passa_alta(fig_folder)
  nivel_mar_filtrado = aplica_filtro(b,a, nivel_mar, tempo, nome_serie, fig_folder)

  return nivel_mar_filtrado

def passa_baixa(fig_folder, modelo):
  # vai ser aplicado em dois passos pra mim: primeior passa baixa c freq horaria
  # e depois passa alta c freq diaria
  if modelo == False:
    dt = 1/24  # Frequência de amostragem (1 hora)
    pesos = 10
  else:
    dt = 1
    pesos = 44
  periodo_min = 3  # Período mínimo em dias
  cutoff_max = 1 / periodo_min
  fn = 1/(2*dt) # frequencia de Nyquist

  b, a = signal.butter(pesos, cutoff_max/fn,
                        btype='low')  # Filtro passa-banda de ordem 2
  frequencias, resposta = signal.freqz(b, a)
  plt.figure()
  plt.semilogx(1/(frequencias*fn/np.pi), np.abs(resposta))
  plt.title('Resposta de Frequência do Filtro Digital Passa Baixa')
  plt.xlabel('Frequência Normalizada')
  plt.ylabel('Magnitude')
  plt.grid(which='both', axis='both')
  plt.savefig(f'{fig_folder}filtro_passa_baixa.png')


  return (b,a)

def passa_alta(fig_folder):
  '''
    Nesse caso, precisa rodar para o ssh já reamostrado em 1 dia
  '''

  dt = 1  # Frequência de amostragem (1 hora)
  periodo = 30  # Período mínimo em dias
  cutoff = 1 / periodo
  fn = 1/(2*dt) # frequencia de Nyquist

  b, a = signal.butter(14, cutoff/fn,
                        btype='high')  # Filtro passa-banda de ordem 2
  frequencias, resposta = signal.freqz(b, a)
  plt.figure()
  plt.semilogx(1/(frequencias*fn/np.pi), np.abs(resposta))
  plt.title('Resposta de Frequência do Filtro Digital Passa Alta')
  plt.xlabel('Frequência Normalizada')
  plt.ylabel('Magnitude')
  # plt.grid()
  plt.grid(which='both', axis='both')
  plt.savefig(f'{fig_folder}filtro_passa_alta.png')


  return (b,a)

def passa_banda(fig_folder, modelo):
  # Filtragem passa-banda entre 3 e 30 dias
  if modelo == False:
    dt = 1/24  # Frequência de amostragem (1 hora)
    pesos = 4
  else:
    dt = 1
    pesos = 15
  periodo_min = 3  # Período mínimo em dias
  periodo_max = 30  # Período máximo em dias
  cutoff_max = 1 / periodo_min
  cutoff_min = 1 / periodo_max
  fn = 1/(2*dt) # frequencia de Nyquist
  b, a = signal.butter(pesos, [cutoff_min/fn, cutoff_max/fn],
                        btype='bandpass')  # Filtro passa-banda de ordem 2
  frequencias, resposta = signal.freqz(b, a)
  plt.figure()
  plt.semilogx(1/(frequencias*fn/np.pi), np.abs(resposta))
  plt.title('Resposta de Frequência do Filtro Digital Passa Banda')
  plt.xlabel('Frequência Normalizada')
  plt.ylabel('Magnitude')
  plt.grid(which='both', axis='both')
  plt.savefig(f'{fig_folder}filtro_passa_banda.png')

  return (b,a)

def aplica_filtro(b,a, nivel_mar, tempo, nome_serie, fig_folder):
  plt.ion()


  # nivel_mar_filtrado = signal.lfilter(b, a, nivel_mar)
  nivel_mar_filtrado = signal.filtfilt(b, a, nivel_mar) # essa funcao 'vai e volta'



# Plot os sinais originais e filtrados
  plt.figure(figsize=(10, 6))
  plt.plot(tempo, nivel_mar, label='Sinal original')
  plt.plot(tempo, nivel_mar_filtrado, label='Sinal filtrado')
  plt.xlabel('Tempo')
  plt.ylabel('Nível do mar')
  plt.title(f'Filtragem de sinal de nível do mar - {nome_serie}')
  plt.legend()
  plt.grid(True)

  plt.savefig(f'{fig_folder}{nome_serie}_filtrado.png')

# Plot os sinais originais e filtrados sobrepostos

  plt.figure(figsize=(10, 6))
  plt.plot(tempo, nivel_mar - nivel_mar.mean(), label='Sinal original - média')
  plt.plot(tempo, nivel_mar_filtrado, label='Sinal filtrado')
  plt.xlabel('Tempo')
  plt.ylabel('Nível do mar')
  plt.title(f'Filtragem de sinal de nível do mar - {nome_serie}')
  plt.legend()
  plt.grid(True)

  plt.savefig(f'{fig_folder}{nome_serie}_filtrado_sobreposto.png')

  plt.close('all')

  return nivel_mar_filtrado

  
