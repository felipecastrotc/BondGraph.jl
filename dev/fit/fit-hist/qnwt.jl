using LinearAlgebra
include("dctm.jl")

function qnwt(f, ∇, x_0; ϵ=1e-10, sup_lim=1e1, k_max=2000)
  #quad Algoritmo de Busca Quadrática
  #   Essse algoritmo de Busca Quadrática foi desenvolvido como parte da disciplina
  #   'Métodos de Otimização de Sistemas' da Faculdade de Engenharia Mecânica
  #   da Universidade Estadual de Campinas.

  #   Autor: Felipe de Castro Teixeira Carvalho
  #   RA: 190823

  # Valores iniciais das variáveis
  k = 1
  erro = Inf
  D = Matrix{Float64}(I, size(x_0,1), size(x_0,1))
  while k <= k_max && erro >= ϵ
    # Determina-se a função data em função de α
    # d = -∇(x), onde d é direção de descida
    # x_aux = x_0 + α*d, obtenção do x auxiliar
    # A função abaixo integra tudo em uma única linha: x_aux = x_0 + α(-∇(x))
    ∇_0 = ∇(x_0)
    d_k = -D*∇_0
    f_k(α) = f(x_0 + α*d_k)
    # Obtêm-se o valor de α pela busca dicotômica
    α = dctm(f_k, [0, sup_lim], ϵ=ϵ)[1]
    # Passo em x
    x_k = x_0 + α*d_k
    # Utilizando o critério de parada erro absoluto relativo
    if norm(x_0) == 0
      erro = norm(∇(x_k))
    else
      erro = norm(∇(x_k))
      # erro = norm(x_k - x_0)/norm(x_0)
      # erro = mean([norm(∇(x_k)) norm(x_k - x_0)/norm(x_0)])
    end
    ∇_k = ∇(x_k)
    s = x_k - x_0
    y = ∇_k - ∇_0
    if (s'*y)[1] <= 0
      D = eye(size(x_0,1))
    else
      # DFP
      D = D + (s*s')/(s'*y) - (D*y*y'*D)/(y'*D*y);
    end
    # Icrementa-se 1 na variável contadora
    k += 1
    # Atualiza o valor de x_0
    x_0 = x_k
  end
  out = x_0;
  return out, k-1, erro
end
