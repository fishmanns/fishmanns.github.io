# Agent Loop

## 一瞥

Agent 反复执行**观察->决策->执行->获取反馈**的执行过程,知道任务完成/触发停止条件/达到资源上限

## suedocode

```
def run_agent(user_message, llm, tools, max_steps=10):
    messages = [
        {"role": "user", "content": user_message}
    ]

    for step in range(max_steps):
        # 1. 请求 LLM 决策
        response = llm(messages, tools=tools)

        # 2. 如果模型已经给出最终答案
        if response.is_final:
            return response.content

        # 3. 执行模型选择的工具
        result = execute_tool(
            response.tool_name,
            response.arguments
        )

        # 4. 把决策和工具结果加入上下文
        messages.append(response.to_message())
        messages.append({
            "role": "tool",
            "name": response.tool_name,
            "content": str(result)
        })

    return "达到最大迭代次数，任务未完成"
```