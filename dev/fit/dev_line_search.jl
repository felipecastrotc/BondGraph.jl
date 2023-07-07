

f(x...) = @. x[1]^2 + sin(x[2])*20 + 10*x[2]
gd(x) = gradient((x) -> f(x), x)[1]
hd(x) = hessian((x) -> f(x), x)




min

x0 = [8.0, 8.0]


d = -hd(x0)*gd(x0)
f̂(α) = f(x0 .+ α*d)

α = 1e-7:0.0000001:1e-1
α = -2:0.01:4
plot(f̂.(10.0.^α), xaxis=:log, yaxis=:log)

dϕ = gd(x0)'*d

ρ = 0.5
c1 = 1e-3
α0 = 1e-7

ϕ0 = f̂(α0)

ϕ < (ϕ0 + c1*α0*dϕ)
α0 = α0*1.1
ϕ = f̂(α0)

f̂(α1)


function dctm(f, range; ϵ=1e-14, l=NaN, k_max=2000)
  #dctm Algoritmo de Busca Dicotômica
  #   Essse algoritmo de Busca Dicotômica foi desenvolvido como parte da disciplina
  #   'Métodos de Otimização de Sistemas' da Faculdade de Engenharia Mecânica
  #   da Universidade Estadual de Campinas.

  #   Autor: Felipe de Castro Teixeira Carvalho
  #   RA: 190823

  # Valores padrão das variáveis
  isnan(l) && (l = 2.1*ϵ)
  # Determina-se um valor inicial do intervalo atual (l_0)
  l_0 = Inf
  k = 0
  while k <= k_max && l_0 > l
      # Calcula-se o ponto médio
      x_M = (range[1] + range[2])/2
      # Calcula-se os valores de λ e ϵ
      λ = x_M - ϵ
      μ = x_M + ϵ
      # Cria-se um novo intervalo baseado em f(λ) e f(μ)
      if f(λ)[1] > f(μ)[1]
          range = [λ range[2]]
      else
          range = [range[1] μ]
      end
      # Adiciona-se 1 na variável contadora
      k += 1
      l_0 = abs(range[2] - range[1])
  end
  erro = l_0
  # Calcula-se a média do intervalo
  out = mean(range)
  return out, k-1, erro
end

function dctmlog(f, range; ϵ=1e-14, l=NaN, k_max=2000)
  #dctm Algoritmo de Busca Dicotômica
  #   Essse algoritmo de Busca Dicotômica foi desenvolvido como parte da disciplina
  #   'Métodos de Otimização de Sistemas' da Faculdade de Engenharia Mecânica
  #   da Universidade Estadual de Campinas.

  #   Autor: Felipe de Castro Teixeira Carvalho
  #   RA: 190823

  # Valores padrão das variáveis
  isnan(l) && (l = 2.1*ϵ)
  # Determina-se um valor inicial do intervalo atual (l_0)
  l_0 = Inf
  k = 0
  while k <= k_max && l_0 > l
      # Calcula-se o ponto médio
      x_M = 10^(sum(log10.(range))/2)
      # Calcula-se os valores de λ e ϵ
      #   λ = x_M - ϵ
      #   μ = x_M + ϵ
      λ = x_M - 10^(log10(x_M)-5)
      μ = x_M + 10^(log10(x_M)-5)
      # Cria-se um novo intervalo baseado em f(λ) e f(μ)
      if f̂(λ)[1] > f̂(μ)[1]
          range = [λ range[2]]
      else
          range = [range[1] μ]
      end
      # Adiciona-se 1 na variável contadora
      k += 1
      l_0 = abs(range[2] - range[1])
  end
  erro = l_0
  # Calcula-se a média do intervalo
  out = mean(range)
  return out, k-1, erro
end

α1, _, _ = dctm(f̂, [1e-7, 4], ϵ=1e-10)

x0 += α1*d