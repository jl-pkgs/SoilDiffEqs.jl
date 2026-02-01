# model_id

- `google/antigravity-gemini-3-flash`
- `google/antigravity-claude-sonnet-4-5`
- `openai/gpt-5.2`
- `minimax-cn-coding-plan/MiniMax-M2.1`
- `zai-coding-plan/glm-4.7`

# 测试方法

`opencode run --model {model_id} "你是什么模型" `，如果能正确返回模型信息就是正常工作的。

```bash
opencode run --model openai/gpt-5.2 "你是什么模型" 
我是 OpenAI 的 `gpt-5.2`（模型标识：`openai/gpt-5.2`）。
```

# GOAL

帮我并行测试这些模型是否正常，记录模型响应时间，最终整理为md表格。

测试失败的模型，循环测试三次，如果均失败，才认定为失败；若有一次成功则跳出测试。
